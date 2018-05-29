__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu"
__project__ = "LLP - DiTaxa"
__website__ = "https://llp.berkeley.edu/ditaxa/"


import sys

sys.path.append('../')
import random
from utility.file_utility import FileUtility
from Bio import SeqIO
from multiprocessing import Pool
import tqdm
from make_representations.cpe_efficient import train_cpe
import sentencepiece as spm
import timeit


class BioCPESegmentTrainMetagenomics:
    '''
        Training the BioCPE segmentation
    '''

    def __init__(self, file_directory, file_extenstion, onlyfiles=[]):
        '''
        :param file_directory: the samples directory
        :param file_extenstion: the file extension fastq or fasta
        :param onlyfiles: filter a list of files
        :param backend: which backend to use
        '''
        print('Segmentation training')
        self.file_directory = file_directory
        self.file_extenstion = file_extenstion
        self.fasta_files, self.filename_mapping = FileUtility.read_fasta_directory(self.file_directory,
                                                                                   self.file_extenstion,
                                                                                   only_files=onlyfiles)
        print(str(len(self.fasta_files)), 'fasta files found in', self.file_directory)

    def generate(self, vocab_size, sample_size, output_dir, num_p=4, backend='Sentencepiece'):
        '''
        :param vocab_size: the size of final vocabulary
        :param sample_size: how many reads from each file
        :param output_dir: where to write the results
        :param num_p: number of cores
        :return:
        '''
        start = timeit.default_timer()
        fasta_files = [(x, sample_size) for x in self.fasta_files]
        corpus = []
        pool = Pool(processes=num_p)
        for ky, v in tqdm.tqdm(pool.imap_unordered(self._get_corpus, fasta_files, chunksize=num_p),
                               total=len(fasta_files)):
            corpus = corpus + v
        pool.close()
        print('Corpus size for training CPE is ', len(corpus))
        if backend == 'Sentencepiece':
            FileUtility.save_list('../tmp/tmp_txt', corpus)
            spm.SentencePieceTrainer.Train(
                '--input=../tmp/tmp_txt --model_prefix=' + output_dir + ' --add_dummy_prefix false --max_sentencepiece_length=512 --model_type=bpe --mining_sentence_size=5000000 --input_sentence_size=10000000 --vocab_size=50000')
            FileUtility.save_list('../tmp/tmp_txt', corpus[0:10])
        elif backend == 'normalbpe':
            train_cpe(corpus, output_dir, vocab_size, output_dir + '_freq')
        print(' The segmentation training took ', timeit.default_timer() - start, ' ms.')

    def _get_corpus(self, file_name_sample, unqiue_reads=True):
        '''
        :param file_name_sample:
        :return:
        '''
        file_name = file_name_sample[0]
        sample_size = file_name_sample[1]
        corpus = []
        if file_name[-1] == 'q':
            for cur_record in SeqIO.parse(file_name, "fastq"):
                corpus.append(str(cur_record.seq).lower())
        else:
            for cur_record in SeqIO.parse(file_name, "fasta"):
                corpus.append(str(cur_record.seq).lower())
        if unqiue_reads:
            corpus = [x.replace(' ', '') for x in corpus]
            corpus = list(set(corpus))
        if sample_size == -1:
            return file_name, corpus
        else:
            return file_name, random.sample(corpus, min(sample_size, len(corpus)))


if __name__ == '__main__':
    G16s = BioCPESegmentTrainMetagenomics('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/dental/',
                                          'fastq')
    G16s.generate(50000, 5000,
                  '/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/dental/test/dental_unique_50000v_5000s',
                  backend='Sentencepiece')
