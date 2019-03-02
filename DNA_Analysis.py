"""
Danny Diaz
eid: dd32387
2-7-19
BCH 394: Bioinformatics

Requires python3.6+ to run
"""



from collections import Counter
from datetime import datetime
from itertools import product
from textwrap import fill
from scipy.stats import chisquare
from contextlib import redirect_stdout

import regex as re




class DNA_Analysis():

    native_nuc_letters = 'ACGT'
    start_codon ='ATG'
    stop_codons = ('TAA', 'TAG', 'TGA')

    def __init__(self,input_DNA, native_nuc_only=True):
        DNA_file = re.compile(r'(/|\\)?.+\.\w{1,10}$')
        self._native_nuc_only = native_nuc_only

        if DNA_file.search(input_DNA):
            self._inp_file = input_DNA
            self.dna_seq = self._get_dna_seq()
        else:
            self.dna_seq = input_DNA

        self._all_nuc_letters = self.get_all_nuc_letters_in_file()

        if self._native_nuc_only:
            self.nucleotides = DNA_Analysis.native_nuc_letters
        else:
            self.nucleotides = self._all_nuc_letters

        self.dna_length = self.calc_genome_length()
        self.mono_nuc_freq = self.set_mono_nuc_freq()
        self.di_nuc_freq = self.set_di_nuc_freq()



    def __repr__(self):
        '''
        String representation on how to recreate the object.
        '''
        return f'DNA_Analysis(file_path={self.file_path}, _native_nuc_only={self.native_nuc_letters})'


    def _get_dna_seq(self):
        with open(self._inp_file, "rt") as inp_file:
            return ''.join(letter for letter in inp_file.read() if not (letter == '\n'))



    def get_gene_region_indices(self,min_peptide_len=50):
        '''
        This will return all possible gene substring indices as a tuple (start_codon, stop_codon).
        To be a match, the gene must transcribe a peptide that is at least 50 AA (default) in length and the start
        and stop codon must be in the same reading frame.
        This is a permutation. It only works on DNA sequences that are of a small size (a few thousands).
        It will take hours on entire genomes.
        '''
        min_nuc_len = min_peptide_len * 30

        gene_start = re.compile('ATG', re.MULTILINE | re.IGNORECASE)
        gene_stop = re.compile('TAA|TAG|TGA', re.MULTILINE | re.IGNORECASE)

        all_start_codons = [codon.start() for codon in gene_start.finditer(self.dna_seq,  overlapped=True)]
        all_stop_codons = [codon.start() for codon in gene_stop.finditer(self.dna_seq,  overlapped=True)]

        return [(start_codon, stop_codon)
                    for start_codon in all_start_codons
                        for stop_codon in all_stop_codons
                            if (stop_codon-start_codon) >= min_nuc_len and (stop_codon-start_codon)%3 == 0]




    def get_all_nuc_letters_in_file(self):
        '''
        Total possible values (characters) for nucleotide positions present in the genome file provided.
        This excludes the newline character (\n) by default.
        '''
        return ''.join({letter for letter in self.dna_seq if not (letter == '\n')})



    def calc_genome_length(self):
        '''
        Total # nucleotides in a Genome. There is the option to include non-native characters into
        the count. The function defaults to only counting native characters (ATCG).
        '''
        if self._native_nuc_only:
            return len([letter for letter in self.dna_seq if letter in 'ATCG'])
        else:
            return len([letter for letter in self.dna_seq if not (letter == '\n')])



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
        if self._native_nuc_only:
            counts = Counter([letter for letter in self.dna_seq if letter in 'ATCG'])
        else:
            counts = Counter([letter for letter in self.dna_seq if not (letter == '\n')])
        return counts


    def set_mono_nuc_freq(self):
        mono_nuc_counts = self.set_mono_nuc_count()
        return {key: round(100 * value / self.dna_length, 3) for key, value in mono_nuc_counts.items()}




    # Question 2 & 3

    def set_di_nuc_count(self,use_native_nuc=False):
        if self._native_nuc_only or use_native_nuc:
            all_dinuc = ["".join(pair) for pair in product(DNA_Analysis.native_nuc_letters, repeat=2)]
        else:
            all_dinuc = ["".join(pair) for pair in product(self._all_nuc_letters, repeat=2)]

        sequence = ''.join(line.strip() for line in self.dna_seq)
        return {dinuc: len(re.findall(dinuc, sequence, overlapped=True))
                       for dinuc in all_dinuc}


    def set_di_nuc_freq(self, use_native_nuc=False):
        '''
        If _native_nuc_only is set to False and use_native_nuc is set to True,
         you will obtain the percent of "other" nucleotides calculated using the total amount of 'ATGC' nucleotides
         present in the genome.
         If _native_nuc_only and use_native_nuc are both set to False then the percent calculated for "other"
         is using the total amount of nucleotides (including other random letters) present in the genome,
         not just the total amount of 'ATGC'.
         It is recommended that use_native_nuc=True unless there is a particular reason. Default assumes you have
         a particular reason.
         '''
        di_nuc_count = self.set_di_nuc_count(use_native_nuc)
        self.di_nuc_freq = {key: round(100 * value / self.dna_length, 3) for key, value in di_nuc_count.items()}
        if use_native_nuc:
            self.di_nuc_freq['other'] = round(100.0 - sum(self.di_nuc_freq.values()), 3)
        return self.di_nuc_freq




    # Question 4

    def exp_di_nuc_freq(self):
        all_dinuc = ["".join(pair) for pair in product(DNA_Analysis.native_nuc_letters, repeat=2)]
        return {dinuc: round(self.ret_exp_freq_di_nuc(dinuc), 3)  for dinuc in all_dinuc}


    # This function does nothing yet.
    def di_nuc_residuals(self):
        exp_freq = self.exp_di_nuc_freq()
        obs_freq = self.di_nuc_freq



    # Question 5

    def is_gene_in_genome(self, gene):
        '''
        Performs a Chi Square Goodness of Fit test on the observed dinucleotide frequencies.
        :param gene:
        :return p-value (str):
        '''
        if isinstance(gene, (str, DNA_Analysis)):
            if isinstance(gene, str):
                gene = DNA_Analysis(gene)
        else:
            raise TypeError('gene argument must be an instance of DNA_Analysis or a '
                            'string of nucleotides.')
        return f'{chisquare(list(gene.di_nuc_freq.values()),list(self.di_nuc_freq.values()))[1]:.2e}'



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



    @classmethod
    def print_HW_1_to_file(cls, file_name='Finished_HW_1.txt'):

        influenzae_native = cls.h_influenzae()
        influenzae_non_native = cls.h_influenzae(native_nuc_only=False)

        aquaticus_native = cls.t_aquaticus()
        aquaticus_non_native = cls.t_aquaticus(native_nuc_only=False)

        mystery_gene_1 = DNA_Analysis('MysteryGene1.txt')
        mystery_gene_2 = DNA_Analysis('MysteryGene2.txt')
        mystery_gene_3 = DNA_Analysis('MysteryGene3.txt')


        with open(file_name,'wt') as f, redirect_stdout(f):

            print('Danny Diaz','eid: dd32387','2-9-19','BCH 394: HW 1',sep='\n')

            print('\n\n', 105 * '*', '\n\n')




            print('Question 1')

            print('\n\n', 105 * '*', '\n\n')

            print('Influenzae Observed Freq:')
            print(fill(f'{influenzae_native.mono_nuc_freq}',108), end='\n\n')

            print('Influenzae Observed Freq: ')
            print(fill(f'{influenzae_non_native.mono_nuc_freq}',100), end='\n')

            print('\n\n', 105 * '*', '\n\n')

            print('Aquaticus Observed Freq:  ')
            print(fill(f'{aquaticus_native.mono_nuc_freq}',108), end='\n\n')

            print('Aquaticus Observed Freq:  ')
            print(fill(f'{aquaticus_non_native.mono_nuc_freq}',108), end='\n')

            print('\n\n', 105 * '*', '\n\n')




            print('Question 2, 3, and 4')

            print('\n\n', 105 * '*', '\n\n')

            print('Influenzae Expected Freq: ')
            print(fill(f'{influenzae_native.exp_di_nuc_freq()}',108), end='\n\n')

            print('Influenzae Observed Freq:   ')
            print(fill(f'{influenzae_native.di_nuc_freq}',108), end='\n')

            print('\n\n',105 * '*', '\n\n')

            print('Aquaticus Expected Freq: ')
            print(fill(f'{aquaticus_native.exp_di_nuc_freq()}',108), end='\n\n')

            print('Aquaticus Observed Freq:   ')
            print(fill(f'{aquaticus_native.di_nuc_freq}',108), end='\n')

            print('\n\n', 105 * '*', '\n\n')

            print(fill('The observed frequencies for both organisms do not match the expected frequencies. I believe '
                       'this is due to the presence of genes in the genome. Genes consist of 3 nucleotide codons that'
                       'are translated into amino acids. The evolutionary pressure on these proteins to conserves their '
                       'sequence is reflected in observed dinuecleotide frequencies observed.', 105))

            print('\n\n',105 * '*', '\n\n')




            print('Question 5')

            print('\n\n', 105 * '*', '\n\n')

            print('Mystery Gene 1 expected Freq: ')
            print(fill(f'{mystery_gene_1.exp_di_nuc_freq()}',108), end='\n\n')

            print('Mystery Gene 1 Observed Freq:   ')
            print(fill(f'{mystery_gene_1.di_nuc_freq}',108), end='\n')

            print('\n\n', 105 * '*', '\n\n')

            print('Mystery Gene 2 expected Freq: ')
            print(fill(f'{mystery_gene_2.exp_di_nuc_freq()}',108), end='\n\n')

            print('Mystery Gene 2 Observed Freq:   ')
            print(fill(f'{mystery_gene_2.di_nuc_freq}',108), end='\n')

            print('\n\n', 105 * '*', '\n\n')

            print('Mystery Gene 3 expected Freq: ')
            print(fill(f'{mystery_gene_3.exp_di_nuc_freq()}',108), end='\n\n')

            print('Mystery Gene 3 Observed Freq:   ')
            print(fill(f'{mystery_gene_3.di_nuc_freq}',108), end='\n')

            print('\n\n', 105 * '*', '\n\n')

            print(f'p-value that mystery gene 1 is in the Influenza genome: '
                  f'{influenzae_native.is_gene_in_genome(mystery_gene_1)}')

            print(f'p-value that mystery gene 2 is in the Influenza genome: '
                  f'{influenzae_native.is_gene_in_genome(mystery_gene_2)}')

            print(f'p-value that mystery gene 3 is in the Influenza genome: '
                  f'{influenzae_native.is_gene_in_genome(mystery_gene_3)}')

            print('\n\n', 105 * '*', '\n\n')

            print(f'p-value that mystery gene 1 is in the Aquaticus genome: '
                  f'{aquaticus_native.is_gene_in_genome(mystery_gene_1)}')

            print(f'p-value that mystery gene 2 is in the Aquaticus genome: '
                  f'{aquaticus_native.is_gene_in_genome(mystery_gene_2)}')


            print(f'p-value that mystery gene 3 is in the Aquaticus genome: '
                  f'{aquaticus_native.is_gene_in_genome(mystery_gene_3)}')

            print('\n\n', 105 * '*', '\n\n')

            print(fill('The p-values provided are from a "Chi-Square Goodness of Fit" test conducted on the observed'
                       'dinucleotide frequencies between each gene and each genome. The comparison of the Hinfluenzae '
                       'genome with mystery gene 1 and 3 provided the respective p-values:.915 and .999. Thus, we cannot '
                       'reject the null hypothesis that their is a statistically significant difference between their '
                       'observed dinucleotide frequencies. Hence, mystery gene 1 and 3 are in the Hinfluenzae genome. '
                       'The p-value for mystery gene 2 was statistically significant (< 0.05). Thus, we reject the null '
                       'hypothesis that there is not a difference between the observed dinucleotide frequencies of '
                       'mystery gene 2 and Hinfluenzae.',
                       105), end='\n\n')

            print(fill('The comparison of the Taquaticus genome with the 3 mystery gene provided 2 p-values that are '
                       'statistically significant (< 0.05): gene 1 and 3. Thus, we can reject the null hypothesis that '
                       'that their is no difference between their observed dinucleotide frequencies. This makes perfect '
                       'sense because we concluded that these genes were present in Hinfluenzae. The p-value for mystery '
                       'gene 2 is 1.000; the null hypothesis cannot be rejected and mystery gene 2 observed dinucleotide '
                       'frequencies show no difference from the the dinucleotide frequencies observed in Taquaticus. '
                       'Hence, mystery gene 2 is in the Taquaticus genome.', 105))

            print('\n\n', 105 * '*', '\n\n')




            print('Question 6')

            print('\n\n', 105 * '*', '\n\n')

            print(fill('You can store your a copy of the entire human genome on your cell phone\'s 100 GB SD card if you '
                       'assume each nucleotide is 1 byte. The human genome is 3 giga base pairs and you have 100 giga '
                       'bytes available.', 105), end='\n\n')

            print(fill('To sequence the genome of every human diagnosed with cancer in the US each year (1,735,350 in '
                       '2018) will cost a little over $1.735 billion, assuming each genome is $1000 and we do not receive '
                       'whole sale prices.', 105), end='\n\n')

            print(fill('To store the genomes of every American diagnosed with cancer in 2018 you will require 5.21 Peta '
                       'bytes of space on your local PC.', 105))

            print('\n\n', 105 * '*', '\n\n')




            print('Question 7')

            print('\n\n', 105 * '*', '\n\n')

            print(fill('The human genome is approximately 652 times larger than the E. Coli genome. ',105), end='\n\n')

            print(fill('The density of genes in the human genome is approximately 133,333 bp per gene. This approximatation'
                       'assumes that there is 22,500 genes in the human genome and that the genome is 3 giga base pairs.'
                       'The density of genes in E. Coli is approximately 1,022 bp per gene. This is assuming that there '
                       'are 4,500 genes and the genome length is 4.6 mega base pairs.', 105))

            print('\n\n', 105 * '*', '\n\n')



            print('Question 8')

            print('\n\n', 105 * '*', '\n\n')

            print(fill('Tryptophan is the amino acid that is least likely to be substituted by another. It has a '
                       'negative substitution score with every amino acid besides phenylalanine and tyrosine, '
                       'the other 2 aromatic amino acids. These scores are 1 and 2 respectively, however, with '
                       'itself its blossum score is 15. Thus, it really does not like to be substituted by other '
                       'aromatic amino acids when compared to self.', 105), end='\n\n')

            print(fill('Other amino acids that show a similar profile to tryptophan are glycine and cysteine. For '
                       'cysteine, its score with itself is 13 compared to tryptophan\'s 15. However, unlike '
                       'tryptophan it does not have a positive blossum score for any other amino acid, but its '
                       'score values are not as negative as tryptophan. Like cysteine, glycine does not have a '
                       'positive blossum score for any substitution except for itself. However, the blossum score '
                       'with itself is only an 8 and its most positive blossum score is 3 zero values with alanine, '
                       'serine, and asparagine.', 105))

            print('\n\n', 105 * '*', '\n\n')



            print('Question 9')

            print('\n\n', 105 * '*', '\n\n')

            print(fill('Polar amino acids are most easily substituted by other polar amino acids. Specifically, '
                       'asparagine and glutamine. Especially if you consider substitution scores of only non-negative '
                       'numbers (including 0), these 2 amino acids scored the highest each with 9. In close second, '
                       'lysine, histidine, and glutamic acid each scored 7 and serine scored 8. ', 105))

            print('\n\n', 105 * '*', '\n\n')



            print('Question 10')

            print('\n\n', 105 * '*', '\n\n')

            print(fill('The most disfavored substitutions are: tryptophan to aspartic acid, tryptophan to cysteine, '
                       'and phenylalanine to aspartic acid. The trend here is the substitution of an aromatic for a '
                       'short acidic residue is highly disfavored.', 105))

            print('\n\n', 105 * '*', '\n\n')



            print('Question 11')

            print('\n\n', 105 * '*', '\n\n')

            print(fill('The average blossum score between DEH: 1/3', 105), end='\n\n')

            print(fill('The average blossum score between VIL: 7/3', 105), end='\n\n')

            print(fill('The average blossum score between the 2 groups: -11/3', 105), end='\n\n')

            print(fill('Substitution between polar amino acids and greasy amino acids is favorable. '
                       'Subsitution from a polar to a nonpolar or vice-versa is not favorable. '
                       'This is due to their chemical properties and roles within a protein. ', 105))

            print('\n\n', 105 * '*', '\n\n')






if __name__=="__main__":
    DNA_Analysis.print_HW_1_to_file()

