__author__ = "Ehsaneddin Asgari"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu or ehsaneddin.asgari@helmholtz-hzi.de"
__project__ = "LLP - Life Language Processing"
__website__ = "https://llp.berkeley.edu/"

import sys

sys.path.append('../')
from sklearn.utils import resample
from sklearn.preprocessing import normalize
from sklearn.feature_extraction.text import TfidfVectorizer
import itertools
import random
import numpy as np
import matplotlib
matplotlib.use('agg')
from utility.math_utility import get_kl_rows
from matplotlib.pyplot import *
from Bio import SeqIO
from utility.file_utility import FileUtility
import matplotlib.pyplot as plt
from make_representations.cpe_apply import CPE
import _pickle as pickle

class BootStrapping(object):
    '''
    Bootstrapping to find a proper sampling size
    '''

    def __init__(self, input_dir, output_dir, sampling_sizes=[10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000],
                 n_resamples=10, seqtype='fasta', M=10):
        '''
        :param input_dir: input directory or a list of files
        :param output_dir: a directory to generate the output files
        :param sampling_sizes: a list of sampling sizes (N's)
        :param n_resamples: number of resamples from each file (N_R = 10)
        :param seqtype: file suffixes fastq or fasta etc.
        :param M: number of files from the directory to make samplings from (M)
        '''
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.M = M
        self.seqtype = seqtype
        self.sampling_sizes = sampling_sizes
        self.n_resamples = n_resamples
        # dictionary for different k-values
        # are the sampling sizes
        self.N_axis = dict()
        # are mean of distances among k-mer distribution resamples
        self.D_S = dict()
        # are std of distances among k-mer distribution resamples
        self.std_D_S = dict()
        # are mean of distances between k-mer distribution resamples and the whole sample k-mer distribution
        self.D_R = dict()
        # are std of distances between k-mer distribution resamples and the whole sample k-mer distribution
        self.std_D_R = dict()

    def add_cpe_sampling(self, cpe_file, cpe_size):
        '''
        :param cpe: cpe bootstrapping to add in the analysis
        :return:
        '''
        self.N_axis[cpe_size], self.D_S[cpe_size], self.std_D_S[cpe_size], self.D_R[cpe_size], self.std_D_R[cpe_size]= self.get_stats_samples_cpe(cpe_file, cpe_size)

    def get_stats_samples_cpe(self, cpe_file, cpe_size):
        '''
        :param cpe_file:
        :param cpe_size:
        :return:
        '''
        x=[]
        y=[]
        y_tot=[]
        error=[]
        error_tot=[]

        f=open(cpe_file,'r')
        CPE_Applier=CPE(f,separator='', merge_size=cpe_size)
        cpe_vocab=[''.join(x.split()).replace('</w>','').lower() for x in FileUtility.load_list(cpe_file)[1::]]
        cpe_vocab=list(set(cpe_vocab))
        cpe_vocab.sort()
        cpe_vectorizer = TfidfVectorizer(use_idf=False, vocabulary=cpe_vocab, analyzer='word',
                                          norm=None, stop_words=[], lowercase=True, binary=False, tokenizer=str.split)
        # To find the files
        if isinstance(self.input_dir, str):
            sample_files=FileUtility.recursive_glob(self.input_dir,"*"+self.seqtype)
        else:
            sample_files=self.input_dir
        sample_files=random.sample(sample_files, self.M)

        # To iterate over the sampling sizes
        for sample_size in self.sampling_sizes:
            distance_i=[]
            tot_dist_i=[]
            print(' sampling size ' ,sample_size , ' is started ...')
            # To iterate over random files
            for sample_file in sample_files:
                comp_dist=self._get_cpe_distribution(sample_file, CPE_Applier, cpe_vectorizer, -1, 1)
                resamples_cpes=self._get_cpe_distribution(sample_file, CPE_Applier, cpe_vectorizer, sample_size, self.n_resamples)
                distance_i.append(np.mean(get_kl_rows(np.array(resamples_cpes))))
                tot_dist_i=tot_dist_i+list(get_kl_rows(np.vstack((np.array(resamples_cpes),comp_dist[0])))[0:10,10])
            print(' sampling size ' ,sample_size ,  ' is completed.')
            mean_distance=np.mean(distance_i)
            std_distance=np.std(distance_i)
            mean_total_distance=np.mean(tot_dist_i)
            std_total_distance=np.std(tot_dist_i)
            x.append(sample_size)
            y.append(mean_distance)
            error.append(std_distance)
            y_tot.append(mean_total_distance)
            error_tot.append(std_total_distance)
        return x, y, error, y_tot, error_tot

    def _get_cpe_distribution(self, sample_file, CPE_Applier, cpe_vectorizer, sample_size, number_of_samples):
        '''
        :param sample_file:
        :param sample_size:
        :param number_of_samples:
        :return:
        '''
        corpus=[]
        for cur_record in SeqIO.parse(sample_file, 'fasta' if self.seqtype=='fsa' or  self.seqtype=='fa' else self.seqtype):
            corpus.append(str(cur_record.seq).lower())
        if sample_size==-1:
            sample_size=len(corpus)
            resamples=[corpus]
        else:
            resamples=[]
            for i in range(number_of_samples):
                resamples.append(resample(corpus, replace=True, n_samples=sample_size))
        vect=[]
        for rs in resamples:
            rs=[CPE_Applier.segment(x) for x in rs]
            vect.append(normalize(np.sum(cpe_vectorizer.fit_transform(rs).toarray(), axis=0).reshape(1, -1), axis=1, norm='l1')[0])
        return vect


    def add_kmer_sampling(self, k_mer):
        '''
        :param k_mer: k_mer bootstrapping to add in the analysis
        :return:
        '''
        self.N_axis[k_mer], self.D_S[k_mer], self.std_D_S[k_mer], self.D_R[k_mer], self.std_D_R[
            k_mer] = self.get_stats_samples(k_mer)

    def get_stats_samples(self, k_mer):
        '''
        get the D_R and D_S
        :param k_mer:
        :return:
        '''
        x = []
        y = []
        y_tot = []
        error = []
        error_tot = []
        # To find the files
        if isinstance(self.input_dir, str):
            sample_files = FileUtility.recursive_glob(self.input_dir, "*" + self.seqtype)
        else:
            sample_files = self.input_dir
        sample_files = random.sample(sample_files, self.M)

        # To iterate over the sampling sizes
        for sample_size in self.sampling_sizes:
            distance_i = []
            tot_dist_i = []
            print(' sampling size ', sample_size, ' is started ...')
            # To iterate over random files
            for sample_file in sample_files:
                comp_dist = self._get_kmer_distribution(sample_file, k_mer, -1, 1)
                resamples_kmers = self._get_kmer_distribution(sample_file, k_mer, sample_size, self.n_resamples)
                distance_i.append(np.mean(get_kl_rows(np.array(resamples_kmers))))
                tot_dist_i = tot_dist_i + list(
                    get_kl_rows(np.vstack((np.array(resamples_kmers), comp_dist[0])))[0:10, 10])
            print(' sampling size ', sample_size, ' is completed.')
            mean_distance = np.mean(distance_i)
            std_distance = np.std(distance_i)
            mean_total_distance = np.mean(tot_dist_i)
            std_total_distance = np.std(tot_dist_i)
            x.append(sample_size)
            y.append(mean_distance)
            error.append(std_distance)
            y_tot.append(mean_total_distance)
            error_tot.append(std_total_distance)
        return x, y, error, y_tot, error_tot

    def _get_kmer_distribution(self, sample_file, k, sample_size, number_of_samples):
        '''
        generate k-mer distribution
        :param sample_file:
        :param k:
        :param sample_size:
        :param number_of_samples:
        :return:
        '''
        vocab = [''.join(xs) for xs in itertools.product('atcg', repeat=k)]
        vectorizer = TfidfVectorizer(use_idf=False, vocabulary=vocab, analyzer='char', ngram_range=(k, k),
                                     norm=None, stop_words=[], lowercase=True, binary=False)
        corpus = []
        for cur_record in SeqIO.parse(sample_file, 'fasta' if self.seqtype == 'fsa' or self.seqtype=='fa' else self.seqtype):
            corpus.append(str(cur_record.seq).lower())
        if sample_size == -1:
            sample_size = len(corpus)
            resamples = [corpus]
        else:
            resamples = []
            for i in range(number_of_samples):
                resamples.append(resample(corpus, replace=True, n_samples=sample_size))
        vect = []
        for rs in resamples:
            vect.append(
                normalize(np.sum(vectorizer.fit_transform(rs).toarray(), axis=0).reshape(1, -1), axis=1, norm='l1')[0])
        return vect

    def plotting(self, file_name, dataset_name):
        '''
        Plotting
        :param file_name:
        :return:
        '''
        figure(figsize=(20, 10))
        matplotlib.rcParams['mathtext.fontset'] = 'stix'
        matplotlib.rcParams['font.family'] = 'STIXGeneral'
        matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
        matplotlib.rcParams['mathtext.fontset'] = 'custom'
        matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
        matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
        matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
        matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
        plt.rc('text', usetex=True)

        ax = subplot(121)
        k_mers = list(self.N_axis.keys())
        legend_vals = []
        for k in k_mers:
            ax.plot(self.N_axis[k], np.array(self.D_S[k]))
            ax.fill_between(self.N_axis[k], np.array(self.D_S[k]) - np.array(self.std_D_S[k]),
                            np.array(self.D_S[k]) + np.array(self.std_D_S[k]), alpha=0.2, linewidth=4,
                            linestyle='dashdot', antialiased=True)
            legend_vals.append('k = ' + str(k))
        ax.legend(legend_vals, loc='upper right', prop={'size': 30}, ncol=1)
        xlabel('Resample size (N)', fontsize=24)
        ylabel(r'$\bar{D_S}(N,k)$', fontsize=24)
        xlim([0, 10000])
        ylim([0, 10])
        title(
            r'(i) \textbf{Self-inconsistency $\bar{D_S}$,} with respect to sample size (N)\\ demonstrated for different k values in the '+dataset_name+' dataset',
            fontsize=24, y=1.01)

        ax = subplot(122)

        k_mers = list(self.N_axis.keys())
        legend_vals = []
        for k in k_mers:
            ax.plot(self.N_axis[k], np.array(self.D_R[k]))
            ax.fill_between(self.N_axis[k], np.array(self.D_R[k]) - np.array(self.std_D_R[k]),
                            np.array(self.D_R[k]) + np.array(self.std_D_R[k]), alpha=0.2, linewidth=4,
                            linestyle='dashdot', antialiased=True)
            legend_vals.append('k = ' + str(k))
        ax.legend(legend_vals, loc='upper right', prop={'size': 30}, ncol=1)
        xlim([0, 10000])
        ylim([0, 0.1])
        xlabel('Resample size (N)', fontsize=24)
        ylabel(r'$\bar{D_R}(N,k)$', fontsize=24)
        title(
            r'(ii) \textbf{Unrepresentativeness $\bar{D_R}$,} with respect to sample size (N)\\ demonstrated for different k values in the '+dataset_name+' dataset',
            fontsize=24, y=1.01)
        plt.tight_layout()
        plt.savefig((self.output_dir + file_name + '.pdf').replace('.pickle',''))
        self.save(file_name.replace('.pickle.pickle','.pickle'))

    def plotting_cpe(self, file_name, dataset_name, ylimv=100, xlimv=20000, ratio=10):
        '''
        Plotting
        :param file_name:
        :return:
        '''
        figure(figsize=(20, 10))
        matplotlib.rcParams['mathtext.fontset'] = 'stix'
        matplotlib.rcParams['font.family'] = 'STIXGeneral'
        matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
        matplotlib.rcParams['mathtext.fontset'] = 'custom'
        matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
        matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
        matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
        matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
        plt.rc('text', usetex=True)

        ax = subplot(121)
        k_mers = list(self.N_axis.keys())
        legend_vals = []
        for k in k_mers:
            ax.plot(self.N_axis[k], np.array(self.D_S[k]))
            ax.fill_between(self.N_axis[k], np.array(self.D_S[k]) - np.array(self.std_D_S[k]),
                            np.array(self.D_S[k]) + np.array(self.std_D_S[k]), alpha=0.2, linewidth=4,
                            linestyle='dashdot', antialiased=True)
            legend_vals.append('\#cpe = ' + str(k))
        ax.legend(legend_vals, loc='upper right', prop={'size': 30}, ncol=1)
        xlabel('Resample size (N)', fontsize=24)
        ylabel(r'$\bar{D_S}(N,\# cpe)$', fontsize=24)
        xlim([0, xlimv])
        ylim([0, ylimv])
        title(
            r'(i) \textbf{Self-inconsistency $\bar{D_S}$,} with respect to sample size (N)\\ demonstrated for different \#cpe in the '+dataset_name+' dataset',
            fontsize=24, y=1.01)

        ax = subplot(122)

        k_mers = list(self.N_axis.keys())
        legend_vals = []
        for k in k_mers:
            ax.plot(self.N_axis[k], np.array(self.D_R[k]))
            ax.fill_between(self.N_axis[k], np.array(self.D_R[k]) - np.array(self.std_D_R[k]),
                            np.array(self.D_R[k]) + np.array(self.std_D_R[k]), alpha=0.2, linewidth=4,
                            linestyle='dashdot', antialiased=True)
            legend_vals.append('\#cpe = ' + str(k))
        ax.legend(legend_vals, loc='upper right', prop={'size': 30}, ncol=1)
        xlim([0, xlimv])
        ylim([0, ylimv/ratio])
        xlabel('Resample size (N)', fontsize=24)
        ylabel(r'$\bar{D_R}(N,\# cpe)$', fontsize=24)
        title(
            r'(ii) \textbf{Unrepresentativeness $\bar{D_R}$,} with respect to sample size (N)\\ demonstrated for different \#cpe in the '+dataset_name+' dataset',
            fontsize=24, y=1.01)
        plt.tight_layout()
        file_name=file_name.replace('.pickle','')
        plt.savefig(self.output_dir + file_name + '.pdf')
        self.save(file_name+'.pickle')


    def load(self,filename):
        f = open(filename, 'rb')
        tmp_dict = pickle.load(f)
        f.close()

        self.__dict__.update(tmp_dict)


    def save(self,filename):
        f = open(self.output_dir + filename, 'wb')
        pickle.dump(self.__dict__, f, 2)
        f.close()


if __name__ == '__main__':
    '''
        test-case
    '''

    BT=BootStrapping('a','b','c')
    BT.load('/mounts/data/proj/asgari/dissertation/git_repos/bio_cpe_data/results/datasetting1/setting1.pickle')
    BT.output_dir='/mounts/data/proj/asgari/dissertation/git_repos/bio_cpe_data/results/datasetting1/'
    BT.plotting_cpe('setting1.pickle','small-effect')
