
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import Counter
from itertools import product
from sklearn.preprocessing import LabelEncoder

import pydot as dot
import numpy as np
import pandas as pd
import regex as re
import networkx as nx



class TM_HMM():

    obs_seq = 'KNSFFFFFFFLIII'

    def __init__(self,
                 sol_seq_file='dna_files/hw2/soluble_sequences.fasta',
                 tm_seq_file='dna_files/hw2/transmembrane_sequences.fasta',
                 state_seq_file='dna_files/hw2/state_sequences.txt'):

        self.sol_seq = sol_seq_file
        self.tm_seq  = tm_seq_file
        self.state_seq = state_seq_file

        self.sol_aa_count, self.sol_seq_len = self.count_AA(self.sol_seq)
        self.tm_aa_count, self.tm_seq_len = self.count_AA(self.tm_seq)

        self.sol_aa_freq = self.calc_AA_em_prob(self.sol_aa_count, self.sol_seq_len)
        self.tm_aa_freq = self.calc_AA_em_prob(self.tm_aa_count, self.tm_seq_len)
        self.create_hmm_graph()





    @property
    def state_seq(self):
        return self._state_seq

    @state_seq.setter
    def state_seq(self, seq_file):
        with open(seq_file, 'r') as f:
            self._state_seq = f.read().strip()
            self._total_states = len(self._state_seq)
            self._calc_state_trans_prob()

    """
       ToDo the values returned from count_AA should be private variables and the freq should be
       calculated directly in the __init__
    """
    def count_AA(self, aa_sequence, ftype='fasta'):
        req = next(SeqIO.parse(aa_sequence, ftype))
        aa_counts = ProteinAnalysis(str(req.seq)).count_amino_acids()
        seq_len = np.sum(list(aa_counts.values()))
        return aa_counts,seq_len



    def calc_AA_em_prob(self, seq_count, seq_len):
        return {key: round(value/seq_len, 5) for key,value in seq_count.items()}


    def _calc_state_init_prob(self):
        self.total_proteins = len(self._state_seq.split('\n'))
        self.state_counts = Counter(state for state in self._state_seq if not state == '\n')
        self.state_init_prob = {
            key:round(value/self.total_proteins,5)
            for key, value in Counter(line[0] for line in self._state_seq.split('\n')).items()
        }



        # self.state_init_prob =  {
        #     key: round(value/self._total_states, 5) for key,value in self.state_counts.items()
        # }


    def _calc_state_digram_count(self):
        all_digram = ["".join(digram) for digram in product('ST', repeat=2)]
        self.digram_counts = {
            digram: len(re.findall(digram,self.state_seq, overlapped=True))
                              for digram in all_digram
        }


    def _calc_state_trans_prob(self):
        try:
            self.state_counts
            self.digram_counts
        except:
            self._calc_state_init_prob()
            self._calc_state_digram_count()
        self.state_trans_prob = {
            key:round(value/self.state_counts['S']
                      if key[0] in 'Ss'
                      else value/self.state_counts['T'], 5)
            for key, value in self.digram_counts.items()
        }


    def _construct_hmm_matrices(self):
        self.hidden_states = ('Soluble', 'Transmembrane')
        self.init_states = pd.Series(self.state_init_prob,
                                     index=self.state_init_prob.keys(),
                                     name='states')

        trans_mat = np.array(list(self.state_trans_prob.values())).reshape(2,2)
        self.trans_df = pd.DataFrame(trans_mat,
                                     columns=self.hidden_states,
                                     index=self.hidden_states)
        # self.transition_df.iloc[0,:] = [value for key,value in self.state_trans_prob.items() if key[0] in 'Ss']
        # self.transition_df.iloc[1,:] = [value for key,value in self.state_trans_prob.items() if key[0] in 'Tt']

        self.emissions = tuple(self.sol_aa_count.keys())

        em_mat = np.array([list(self.sol_aa_freq.values()),
                           list(self.tm_aa_freq.values())])

        self.em_df = pd.DataFrame(em_mat, columns=self.emissions, index=self.hidden_states)


    def _create_graph_edges(self):
        self.state_edges = {
            (ind,col):self.trans_df.loc[ind,col]
                for col in self.trans_df.columns
                    for ind in self.trans_df.index
        }

        self.em_edges = {
            (ind,col):self.em_df.loc[ind,col]
                for col in self.em_df.columns
                    for ind in self.em_df.index
        }

    def create_hmm_graph(self):
        try:
            self.em_df
        except:
            self._construct_hmm_matrices()
            self._create_graph_edges()

        self.G = nx.MultiDiGraph()

        self.G.add_nodes_from(self.hidden_states)

        for edges in (self.state_edges, self.em_edges):
            for (frm, to), prob  in edges.items():
                self.G.add_edge(frm, to, weight=prob, label=prob)

        pos = nx.drawing.nx_pydot.graphviz_layout(self.G,prog='neato')
        nx.draw_networkx(self.G, pos)

        pydot_filename = './python_output/hw2_output/transmembrane_aa_hmm.dot'
        png_filename = './python_output/hw2_output/transmembrane_aa_hmm.png'

        nx.drawing.nx_pydot.write_dot(self.G, pydot_filename)
        (graph,) = dot.graph_from_dot_file(pydot_filename)
        graph.write_png(png_filename)



    def run_viterbi(self,obs_seq=174.3):

        aa_encoder = LabelEncoder()
        aa_labels = aa_encoder.fit_transform(self.em_df.columns)
        aa_num_map = {key: value for key, value in zip(self.em_df.columns, aa_labels)}
        num_aa_map = {key: value for key, value in zip(aa_labels, self.em_df.columns)}

        state_encoder = LabelEncoder()
        state_labels = state_encoder.fit_transform(self.init_states.index)
        num_state_map = {key:value for key,value in zip(state_labels, self.init_states.index)}

        if obs_seq == 174.3:
            seq = TM_HMM.obs_seq
        else:
            seq = obs_seq

        seq_map = np.array([aa_num_map[aa] for aa in seq])
        seq_len = len(seq)

        self.init_prob = self.init_states.values
        self.trans_prob = self.trans_df.values
        self.em_prob = self.em_df.values

        nStates = len(self.hidden_states)

        self.path = np.chararray(seq_len)
        ipath = np.zeros(seq_len,dtype=int)

        delta = np.zeros((nStates,seq_len))
        viterbi = np.zeros((nStates,seq_len))

        print('\n\nPrinting all values being used in viterbi:')
        print('State map:', num_state_map, sep='\n', end='\n')
        print('Observed seq:', seq, sep='\n', end='\n')
        print('Obs_seq_map', seq_map, sep='\n', end='\n')
        print('Initial prob:', self.init_states, sep='\n', end='\n')
        print('Transition prob:', self.trans_df, sep='\n', end='\n')
        print('Emission probability', self.em_df.head(), sep='\n', end='\n')
        print('Emission prob of K (initial aa):', self.em_df['K'], sep='\n', end='\n')





        delta[:,0] = self.init_prob * self.em_df.loc[:,seq[0]]
        viterbi[:,0] = 0

        for pos in range(1, seq_len):
            for state in range(nStates):
                print('delta:', delta[:, pos - 1])
                print('trans_prob:', self.trans_prob[:, state])
                print('product: ', delta[:, pos-1] * self.trans_prob[:, state])
                print('em_prob:', self.em_prob[state, seq_map[pos]])
                delta[state,pos] = np.max(
                    delta[:, pos-1] * self.trans_prob[:, state]
                ) * self.em_prob[state, seq_map[pos]]
                print('Max(*em_prob):', delta[state,pos])
                viterbi[state, pos] = np.argmax(
                    delta[:, pos-1] * self.trans_prob[:, state]
                )
                print('max_position:', viterbi[state, pos])

                print(f'state={state}, pos={pos}: viterbi[{state},{pos}] = {viterbi[state,pos]}')

        print(viterbi)


        print('Starting Backtrace')
        #print(seq_len)
        ipath[seq_len-1] = np.argmax(delta[:,seq_len-1])
        #print(ipath)

        for pos in range(seq_len-2, -1, -1):
            ipath[pos] = viterbi[ipath[pos+1],[pos+1]]
            #print(f'ipath[{pos}] = {ipath[pos]}')

        path = ''.join(num_state_map[pos] for pos in ipath)
        print(path)


    def _calc_avg_tm_segments(self):
        num_tm_regions = [len(re.findall('S[T]+?S',line)) for line in self._state_seq.split('\n')]
        avg_tm_region = np.mean(num_tm_regions)
        print(avg_tm_region)



if __name__ == "__main__":
    tm_obj = TM_HMM()
    print(tm_obj.sol_aa_freq)
    print(tm_obj.tm_aa_freq)
    print(tm_obj.sol_aa_count)
    print(tm_obj.state_init_prob)
    print(tm_obj.state_trans_prob.values())
    print(tm_obj.init_states)
    print(tm_obj.trans_df)
    print(tm_obj.em_df)
    print(tm_obj.run_viterbi())
    print(tm_obj._calc_avg_tm_segments())

#
#
# tm_obj.gen_em_prob()
# print(tm_obj.sol_seq.mono_nuc_freq)
# print(tm_obj.tm_seq.mono_nuc_freq)



