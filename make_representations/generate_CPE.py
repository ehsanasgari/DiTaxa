__author__ = "Ehsaneddin Asgari"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu or ehsaneddin.asgari@helmholtz-hzi.de"
__project__ = "LLP - BioCPE"
__website__ = "https://llp.berkeley.edu/micropheno/"

import sys

sys.path.append('../')
import random
from utility.file_utility import FileUtility
from Bio import SeqIO
from multiprocessing import Pool
import tqdm
from make_representations.cpe_efficient import train_cpe


class Generate16SCPE:
    def __init__(self, file_directory, file_extenstion, onlyfiles=[]):
        print('GENERATE CPE')
        self.file_directory = file_directory
        self.file_extenstion = file_extenstion
        self.fasta_files, self.filename_mapping = FileUtility.read_fasta_directory(self.file_directory,
                                                                                   self.file_extenstion,
                                                                                   only_files=onlyfiles)
        print(str(len(self.fasta_files)), 'fasta files found in', self.file_directory)

    def generate(self, vocab_size, sample_size, output_dir, num_p=4):
        fasta_files = [(x, sample_size) for x in self.fasta_files]
        corpus = []
        pool = Pool(processes=num_p)
        for ky, v in tqdm.tqdm(pool.imap_unordered(self.get_corpus, fasta_files, chunksize=5),
                               total=len(fasta_files)):
            corpus = corpus + v
        print('Corpus size for training CPE is ', len(corpus))
        train_cpe(corpus, output_dir, vocab_size, output_dir + '_freq')

    def get_corpus(self, file_name_sample):
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
        if sample_size == -1:
            return file_name, corpus
        else:
            return file_name, random.sample(corpus, min(sample_size, len(corpus)))


if __name__ == '__main__':
    G16s = Generate16SCPE('/mounts/data/proj/asgari/dissertation/git_repos/bio_cpe_data/synthetic/datasetting2', 'fa')
    G16s.generate(50000, 10000,
                  '/mounts/data/proj/asgari/dissertation/git_repos/bio_cpe_data/results/synthetic_setting2_cpe_10000perfile')
