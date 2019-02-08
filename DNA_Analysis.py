from collections import Counter
from datetime import datetime
from itertools import product
from difflib import SequenceMatcher
import regex as re


class DNA_Analysis():

    native_nuc_letters = 'ACGT'

    def __init__(self,file_path, native_nuc_only=True):
        self.inp_file = file_path
        self.native_nuc_only = native_nuc_only
        self.all_nuc_letters = self.get_all_nuc_letters_in_file()
        self.dna_length = self.calc_genome_length()
        self.mono_nuc_freq = self.set_mono_nuc_freq()
        self.di_nuc_freq = self.set_di_nuc_freq()



    def __repr__(self):
        '''
        String representation on how to recreate the object.
        '''
        return f'DNA_Analysis(file_path={self.file_path}, native_nuc_only={self.native_nuc_letters})'


    def get_dna_seq(self):
        with open(self.inp_file, "rt") as inp_file:
            if self.native_nuc_only:
                return ''.join(letter for letter in inp_file.read() if letter in 'ATCG')
            else:
                return ''.join(letter for letter in inp_file.read() if not (letter == '\n'))

    def get_all_nuc_letters_in_file(self):
        '''
        Total possible values (characters) for nucleotide positions present in the genome file provided.
        This excludes the newline character (\n) by default.
        '''
        with open(self.inp_file, "rt") as inp_file:
            return ''.join({letter for letter in inp_file.read() if not (letter == '\n')})



    def calc_genome_length(self):
        '''
        Total # nucleotides in a Genome. There is the option to include non-native characters into
        the count. The function defaults to only counting native characters (ATCG).
        '''
        with open(self.inp_file, "rt") as inp_file:
            if self.native_nuc_only:
                return len([letter for letter in inp_file.read() if letter in 'ATCG'])
            else:
                return len([letter for letter in inp_file.read() if not (letter == '\n')])



    def ret_exp_freq_di_nuc(self,dinuc):
        '''
        Helper function for expected dinucleotide frequency calculations.
        This function is utilizing the mono_nuc_freq from this instance to calculate the expected dinucleotide
        frequency, assuming mutual exclusivity between adjacent nucleotides.
        The probability must be expressed in terms of percentage (0-100%) and not from 0-1 for this function
        to work properly.
        The dinucleotide parameter must be a STRING of 2 nucleotides(characters).
        '''
        return self.mono_nuc_freq[dinuc[0]]*self.mono_nuc_freq[dinuc[1]]/100


    # Question 1

    def set_mono_nuc_count(self):
        with open(self.inp_file, "rt") as inp_file:
            if self.native_nuc_only:
                counts = Counter([letter for letter in inp_file.read() if letter in 'ATCG'])
            else:
                counts = Counter([letter for letter in inp_file.read() if not (letter == '\n')])
            return counts


    def set_mono_nuc_freq(self):
        mono_nuc_counts = self.set_mono_nuc_count()
        return {key: round(100 * value / self.dna_length, 3) for key, value in mono_nuc_counts.items()}




    # Question 2 & 3

    def set_di_nuc_count(self,use_native_nuc=False):
        if self.native_nuc_only or use_native_nuc:
            all_dinuc = ["".join(pair) for pair in product(DNA_Analysis.native_nuc_letters, repeat=2)]
        else:
            all_dinuc = ["".join(pair) for pair in product(self.all_nuc_letters, repeat=2)]
        with open(self.inp_file, 'rt') as inp_file:
            sequence = ''.join(line.strip() for line in inp_file.read())
            return {dinuc: len(re.findall(dinuc, sequence, overlapped=True))
                           for dinuc in all_dinuc}


    def set_di_nuc_freq(self, use_native_nuc=False):
        '''
        If native_nuc_only is set to False and use_native_nuc is set to True,
         you will obtain the percent of "other" nucleotides calculated using the total amount of 'ATGC' nucleotides
         present in the genome.
         If native_nuc_only and use_native_nuc are both set to False then the percent calculated for "other"
         is using the total amount of nucleotides (including other random letters) present in the genome,
         not just the total amount of 'ATGC'.
         It is recommended that use_native_nuc=True unless there is a particular reason. Default assumes you have
         a particular reason.
         '''
        di_nuc_count = self.set_di_nuc_count(use_native_nuc)
        di_nuc_freq = {key: round(100 * value / self.dna_length, 3)
                      for key, value in di_nuc_count.items()}
        if use_native_nuc:
            di_nuc_freq['other'] = round(100.0 - sum([value for value in di_nuc_freq.values()]), 3)
        else:
            self.di_nuc_freq = di_nuc_freq
        return di_nuc_freq




    # Question 4

    def exp_di_nuc_freq(self):
        all_dinuc = ["".join(pair) for pair in product(DNA_Analysis.native_nuc_letters, repeat=2)]
        return {dinuc: round(self.ret_exp_freq_di_nuc(dinuc), 3)  for dinuc in all_dinuc}



    # Question 5
    #ToDo I need to figure out why the sequence matcher is not matching any substrings.
    def gene_in_dna(self,gene):
        if isinstance(gene, DNA_Analysis):
            match = SequenceMatcher(a=gene.get_dna_seq(), b=self.get_dna_seq())
        elif isinstance(gene, str):
            match = SequenceMatcher(a=gene, b=self.get_dna_seq())
        else:
            raise TypeError('gene argument must be an instance of DNA_Analysis or a '
                            'string of nucleotides.')
        longest_match = match.get_matching_blocks()

        print('Longest match length:',match.ratio())
        #print('Longest match gene %:', len(longest_match)/gene.dna_length)



    @classmethod
    def h_influenzae(cls, file_path='Hinfluenzae.txt', native_nuc_only=True):
        return cls(file_path, native_nuc_only)



    @classmethod
    def t_aquaticus(cls,file_path='Taquaticus.txt', native_nuc_only=True):
        return cls(file_path, native_nuc_only)




    #ToDo turn this into a static method or refactor it to spit out only the data of 1 genome at a time.
    def to_file(self,output_filename='GenomeAnalysis_Output'):
        now = datetime.today().strftime('%-m_%-d_%y.%-H:%-M')
        with open(f'{output_filename}_{now}.txt', 'a+') as out_file:
            # Question 1
            flu_data = self.mono_nuc_freq
            aquaticus_data = self.mono_nuc_freq('Taquaticus.txt')
            out_file.write(f'Question 1:\n\tH.Influenzae frequencies: {flu_data}'
                           f'\n\tT.Aquaticus frequencies:  {aquaticus_data}\n\n\n')

            # Question 2
            flu_data = self.di_nuc_freq
            aquaticus_data = self.di_nuc_freq('Taquaticus.txt')
            out_file.write(f'Question 2:\n\tH.Influenzae frequencies: {flu_data}'
                           f'\n\tT.Aquaticus frequencies:  {aquaticus_data}\n\n\n')


if __name__=='__main__':
    # influenzae_native = DNA_Analysis('Hinfluenzae.txt')
    # influenzae_non_native = DNA_Analysis('Hinfluenzae.txt', native_nuc_only=False)

    influenzae_native = DNA_Analysis.h_influenzae()
    # influenzae_non_native = DNA_Analysis.h_influenzae(native_nuc_only=False)

    aquaticus_native = DNA_Analysis.t_aquaticus()
    # aquaticus_non_native = DNA_Analysis.t_aquaticus(native_nuc_only=False)

    mystery_gene_1 = DNA_Analysis('MysteryGene1.txt')

    # This is for question 4 & 5
    print('Influenzae Expected Freq:', influenzae_native.exp_di_nuc_freq())
    print('Influenzae Actual Freq:  ', influenzae_native.di_nuc_freq)

    print('\n\n',240*'*', '\n\n')

    print('Aquaticus Expected Freq:', aquaticus_native.exp_di_nuc_freq())
    print('Aquaticus Actual Freq:  ', aquaticus_native.di_nuc_freq)

    print('\n\n',240*'*', '\n\n')

    print("Mystery Gene 1 expected Freq:", mystery_gene_1.exp_di_nuc_freq())
    print("Mystery Gene 1 Actual Freq:  ", mystery_gene_1.di_nuc_freq)

    print('\n\n', 240 * '*', '\n\n')

    print(influenzae_native.dna_length)
    print(mystery_gene_1.dna_length)
    influenzae_native.gene_in_dna(mystery_gene_1)

    # print(influenzae_native.mono_nuc_freq)
    # print(influenzae_non_native.mono_nuc_freq)
    # print(influenzae_native.di_nuc_freq)
    # print(influenzae_non_native.di_nuc_freq)
    #
    # print('\n\n',180*'*', '\n\n')
    #
    #
    # print(aquaticus_native.mono_nuc_freq)
    # print(aquaticus_non_native.mono_nuc_freq)
    # print(aquaticus_native.di_nuc_freq)
    # print(aquaticus_non_native.di_nuc_freq)