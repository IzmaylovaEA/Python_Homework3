
# coding: utf-8

# In[1]:

from Bio import SeqIO
from collections import defaultdict
import itertools
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import splrep, splev
class Kmer_spectrum:
    def __init__(self, fastq):
        self.fastq = fastq
        self.kmer_dict = defaultdict(int)
    # Извлечение из FASTQ-файла k-меров заданной длины и качества.
    def kmers(self, k, q):
        for record in SeqIO.parse(self.fastq, 'fastq'):
            seq_lng = len(record.seq)
            for index in range(seq_lng-k+1):
                if all(qual >= q for qual in record.letter_annotations['phred_quality'][index:(index+k)]):
                    current_kmer = record.seq[index:(index+k)]
                    current_kmer_rev_compl = current_kmer.reverse_complement()
                    self.kmer_dict[str(current_kmer)] += 1
                    self.kmer_dict[str(current_kmer_rev_compl)] += 1
        return self.kmer_dict
    # Преобразование словаря k-меров в массив.
    def kmer_freq(self, inp_dict):
        return [(x, len(list(y))) for x, y in itertools.groupby(sorted(inp_dict.values()))]
    # Визуализация спектра k-меров.
    def visualize_spectrum(self, inp_list):
        plt.figure(figsize=(50, 25))
        for coordinates in inp_list:
            x = (coordinates[0], coordinates[0])
            y = (0, coordinates[1])
            plt.plot(x, y, color = 'darkblue', marker = 'o')
        inp_list_array = np.array(inp_list)
        x_line = inp_list_array[:,0]
        y_line = inp_list_array[:,1]
        # Сглаживание графика функции.
        bspl = splrep(x_line, y_line, s = 300000)
        y_new = splev(x_line, bspl)
        plt.plot(x_line, y_new, linewidth = 5, color = 'navy')
        # Поиск минимума (граница шумовой и информационной частей спектра) и максимума (основной пик в спектре).
        first_max = y_new[(np.diff(np.sign(np.diff(y_new))) < 0).nonzero()[0]+1][0]
        first_min = y_new[(np.diff(np.sign(np.diff(y_new))) > 0).nonzero()[0]+1][0]
        pos_max = (np.abs(y_line-first_max)).argmin()
        pos_min = (np.abs(y_line-first_min)).argmin()
        plt.axvline(x = x_line[pos_min], linewidth = 7, color = 'darkred', label = 'cut-off')
        plt.legend(fontsize = 50)
        # Преобразование масштаба по оси y в логарифмический вид для более удобного отображения.
        plt.yscale('log')
        plt.xlabel('K-mer multiplicity', fontsize = 50)
        plt.xticks(fontsize= 30)
        plt.xlim(xmin = 0)
        plt.ylabel('Number of distinct k-mers, lg', fontsize = 50)
        plt.yticks(fontsize= 30)
        plt.title('K-mer spectrum', fontsize = 60, fontweight = 'bold')
        plt.show()
        return x_line, y_line, pos_min, pos_max
    # Оценка размера генома.
    def genome_size(self, x_line, y_line, pos_min, pos_max):
        s = 0
        for index in range(pos_min+1, len(x_line)):
            s += x_line[index]*y_line[index]
        size = s/x_line[pos_max]
        print('Genome size is', size, 'bp')
my_file = Kmer_spectrum('/home/izmaylova/Downloads/test_kmer1.fastq')
frequency = my_file.kmer_freq(my_file.kmers(15, 30))
inp = my_file.visualize_spectrum(frequency)
x_line = inp[0]
y_line = inp[1]
pos_min = inp[2]
pos_max = inp[3]
my_file.genome_size(x_line, y_line, pos_min, pos_max)


# In[ ]:



