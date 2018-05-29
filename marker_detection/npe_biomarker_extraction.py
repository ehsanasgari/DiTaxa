__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu"
__project__ = "LLP - DiTaxa"
__website__ = "https://llp.berkeley.edu/ditaxa/"

import sys
sys.path.append('../')
from utility.file_utility import FileUtility
from utility.list_set_util import get_max_of_dict
import numpy as np
from feature_selection.chi2analysis import Chi2Analysis
from scipy.sparse import csr_matrix

class NPEMarkerDetection(object):
    '''
    Class to extract markers to fasta files
    '''
    def __init__(self, X_file, Y_file, features_file, path, selected_samples):
        '''
        :param X:
        :param Y:
        :param features:
        :param path:
        '''
        self.X=FileUtility.load_sparse_csr(X_file)
        self.X=self.X.toarray()
        self.X=self.X[selected_samples,:]
        self.X=csr_matrix(self.X)
        self.Y=[int(x) for x in FileUtility.load_list(Y_file)]
        self.features=FileUtility.load_list(features_file)
        self.path=path


    def extract_markers(self, topn=50000, threshold=5):
        '''
        :param path:
        :param X:
        :param Y:
        :param features:
        :return:
        '''
        CA=Chi2Analysis(self.X, self.Y, self.features)
        res_bin=CA.extract_features_fdr(self.path+'_chi2_binary.txt',topn,direction=True,binarization=True, allow_subseq=False)
        res_med=CA.extract_features_fdr(self.path+'_chi2_median.txt',topn,direction=True,binarization='median', allow_subseq=False)
        pos_bin,neg_bin=self.extract_top_but_n(res_bin,topn,threshold)
        pos_med,neg_med=self.extract_top_but_n(res_med,topn,threshold)
        NPEMarkerDetection.write_in_fastafile(self.path+'_chi2_binary.fasta',res_bin)
        NPEMarkerDetection.write_in_fastafile(self.path+'_chi2_relative.fasta',res_med)
        NPEMarkerDetection.write_in_file(self.path+'_chi2_binary_taxa.txt',pos_bin,neg_bin)
        NPEMarkerDetection.write_in_file(self.path+'_chi2_median_taxa.txt',pos_med,neg_med)

    @staticmethod
    def write_in_fastafile(filename,res, min_length=50):
        corpus=[]
        labels=[]
        for seq, score, pval, _, _ in res:
            if len(seq)>min_length and pval<0.05:
                corpus.append(seq)
                labels.append(' '.join(['+' if score>0 else '-','p-val:'+str(pval)]))
        FileUtility.create_fasta_file(filename,corpus,labels)

    @staticmethod
    def write_in_file(filename,pos,neg):
        lines=[['direction','marker','p-value']]
        for marker, pval in pos:
            lines.append(['+',marker, str(pval)])
        for marker, pval in neg:
            lines.append(['-',marker, str(pval)])
        FileUtility.save_list(filename, ['\t'.join(line) for line in lines])

    def extract_top_but_n(self, chi2_res, n, threshold=3.8):
        # 5
        pos=[]
        neg=[]
        idx=0
        while len(pos)<n and idx < len(chi2_res):
            if chi2_res[idx][1]>threshold :
                pos.append((chi2_res[idx][0], chi2_res[idx][2]))
            idx+=1
        idx=0
        while len(neg)<n and idx < len(chi2_res):
            if chi2_res[idx][1]<-threshold and idx< len(chi2_res):
                neg.append((chi2_res[idx][0], chi2_res[idx][2]))
            idx+=1
        return pos, neg

if __name__=='__main__':
    G16s = NPEMarkerDetection('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/dental/new_cpe_piece/dental_uniqe_seqpiece.npz','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/dental/new_cpe_piece/dental_label.txt','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/dental/new_cpe_piece/dental_uniqe_seqpiece_features','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/dental/new_cpe_piece/dental')
    G16s.extract_markers()
