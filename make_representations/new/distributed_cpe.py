import sys
sys.path.append('../')
from utility.fasta_reader import FastaReader
from representation_utility.cpe_efficient import train_cpe, train_cpe_one_step
from representation_utility.cpe_apply import BPE
from sklearn.utils import resample
from multiprocessing import Pool
import time


class DistributedCPE:
    def __init__(self, sequences, n_samples, sample_size, vocab_size, proc=30,starting_seed=None):
        # init
        self.sequences=sequences
        self.n_samples=n_samples
        self.sample_size=sample_size
        self.proc=proc
        self.vocab_size=vocab_size
        # starting seed
        self.starting_seed=starting_seed

    def cpe_apply(self, sample):
        pool = Pool(self.proc)
        return list(pool.map(self.cpe_applier.segment, sample))

    def one_iter(self, sample):
        return train_cpe_one_step(sample)

    def preprocess(self):
        #sampling
        if self.starting_seed is not None:
            f=open(self.starting_seed,'r')
            self.cpe_applier=BPE(f,'')
            f.close()
        else:
            self.cpe_applier=None
        self.samples=[]
        for i in range(0,self.n_samples):
            if self.cpe_applier is not None:
                self.samples.append(self.cpe_apply(resample(self.sequences, replace=True, n_samples=self.sample_size)))
            else:
                self.samples.append(resample(self.sequences, replace=True, n_samples=self.sample_size))

    def vote_among_samples(self):
        pool = Pool(self.proc)
        self.candidates=list(pool.map(self.one_iter, self.samples))
        with open(self.starting_seed, "a+") as myfile:
            myfile.write(DistributedCPE.most_common(self.candidates))

    def run(self):
        for i in range(self.vocab_size):
            start_time = time.time()
            self.preprocess()
            self.vote_among_samples()
            if i%100==0:
                print ('vocab_size ', str(i))
                print(time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
                start_time = time.time()
    @staticmethod
    def most_common(lst):
        return max(set(lst), key=lst.count)

if __name__ == '__main__':

    FR=FastaReader('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/SILVA/SILVA_128_SSURef_tax_silva.fasta')
    DCPE=DistributedCPE(FR.seq, 10, 100, 100000,50, '/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/16s_start_seed.txt')
    DCPE.run()


