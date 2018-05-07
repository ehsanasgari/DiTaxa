__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu"
__project__ = "LLP - BioCPE"
__website__ = "https://llp.berkeley.edu/16scpe/"

import sys
sys.path.append('../')
from sklearn.preprocessing import normalize
from sklearn.feature_extraction.text import TfidfVectorizer
import itertools
import numpy as np
from multiprocessing import Pool
import tqdm
import random
from scipy import sparse
from utility.file_utility import FileUtility
from Bio import SeqIO
import timeit
from make_representations.cpe_apply import CPE
import sentencepiece as spm

class BioCPESegmentApplyMetagenomics:
    '''
        Make k-mer from directory of fasta files
    '''

    def __init__(self, file_directory, file_extenstion, cpe_file, onlyfiles=[], sampling_number=3000, num_p=20, vocab_size=-1):
        '''
        :param fasta_files: list of fasta files
        :param indexing: the index
        :param sampling_number:
        :param num_p:
        '''
        self.file_directory = file_directory
        self.file_extenstion = file_extenstion
        self.fasta_files, self.indexing = FileUtility.read_fasta_directory(self.file_directory,
                                                                                   self.file_extenstion,
                                                                                   only_files=onlyfiles)
        print(str(len(self.fasta_files)), 'fasta files found in', self.file_directory)

        self.num_p=num_p
        self.sampling_number=sampling_number
        self.cpe_file=cpe_file

        if '.model' in cpe_file:
            self.model_type='seqpiece'
            self.cpe_vocab=[x.split()[0] for x in FileUtility.load_list(cpe_file.replace('.model','.vocab'))]
        else:
            self.model_type='normal_bpe'
            self.cpe_vocab=[''.join(x.split()).replace('</w>','').lower() for x in FileUtility.load_list(cpe_file)[1::]]
            self.cpe_vocab=list(set(self.cpe_vocab))
        self.vocab_size=vocab_size
        self.cpe_vocab.sort()
        self.cpe_vectorizer = TfidfVectorizer(use_idf=False, vocabulary=self.cpe_vocab, analyzer='word',
                                          norm=None, stop_words=[], lowercase=True, binary=False, tokenizer=str.split)

    def generate_cpes_all(self, save=False, norm=False):
        data = np.zeros((len(self.fasta_files), len(self.cpe_vocab))).astype(np.float64)

        # multi processing extraction of cpe distributions
        t_steps=[]
        s_steps=[]
        pool = Pool(processes=self.num_p)
        for ky, (v,t,s) in tqdm.tqdm(pool.imap_unordered(self._get_cpe_distribution, self.fasta_files, chunksize=self.num_p),
                               total=len(self.fasta_files)):
            data[self.indexing[ky], :] = v
            t_steps.append(t)
            s_steps.append(s)
        pool.close()
        # normalize the frequencies
        if norm:
            data = normalize(data, axis=1, norm='l1')
        data = sparse.csr_matrix(data)

        if save:
            FileUtility.save_sparse_csr(save, data)
            FileUtility.save_list(save+'_meta',self.fasta_files)
            FileUtility.save_list(save+'_features',self.cpe_vocab)
            FileUtility.save_list(save+'_log',[': '.join(['mean_time', str(np.mean(t_steps))]), ': '.join(['std_time', str(np.std(t_steps))]), ': '.join(['mean_size', str(np.mean(s_steps))]), ': '.join(['std_size', str(np.std(s_steps))])])
        return data


    def _get_cpe_distribution(self, file_name, make_unique=True):
        '''
        calling it for a single class
        :param file_name:
        :param make_unique:
        :return:
        '''
        if self.model_type=='seqpiece':
            self.CPE_Applier = spm.SentencePieceProcessor()
            self.CPE_Applier.Load(self.cpe_file)
        else:
            f=open(self.cpe_file,'r')
            self.CPE_Applier=CPE(f,separator='', merge_size=self.vocab_size)


        start = timeit.timeit()
        corpus=[]
        if file_name[-1]=='q':
            for cur_record in SeqIO.parse(file_name, "fastq"):
                corpus.append(str(cur_record.seq).lower())
        else:
            for cur_record in SeqIO.parse(file_name, "fasta"):
                corpus.append(str(cur_record.seq).lower())
        if make_unique:
            corpus=list(set(corpus))
        tot_size=len(corpus)
        if self.sampling_number==-1:
            random.shuffle(corpus)
        else:
            corpus = random.sample(corpus, min(self.sampling_number,len(corpus)))
        corpus=[' '.join(self.CPE_Applier.EncodeAsPieces(x)) for x in corpus]
        end = timeit.timeit()
        return file_name,(np.sum(self.cpe_vectorizer.fit_transform(corpus).toarray(), axis=0),end - start,tot_size)


if __name__=='__main__':
    G16s = BioCPESegmentApplyMetagenomics('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/dental/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/dental/new_cpe_piece/dental_unique_50000v_5000s.model',sampling_number=-1)
    G16s.generate_cpes_all('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/dental/new_cpe_piece/dental_uniqe_seqpiece')

