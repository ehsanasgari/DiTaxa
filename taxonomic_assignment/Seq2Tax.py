__author__ = "Ehsaneddin Asgari"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu or ehsaneddin.asgari@helmholtz-hzi.de"
__project__ = "LLP - Life Language Processing"
__website__ = "https://llp.berkeley.edu/"

import sys
sys.path.append('../')
from utility.file_utility import FileUtility
from utility.list_set_util import get_max_of_dict
import numpy as np
from chi2analysis.chi2analysis import Chi2Analysis


class Seq2Tax(object):

    def __init__(self):
        self.seq2tax=FileUtility.load_obj('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/GG/to_be_checked/seq2tax.pickle')


    def extract_markers(self, path, X, Y, features):
        '''
        :param path:
        :param X:
        :param Y:
        :param features:
        :return:
        '''
        CA=Chi2Analysis(X, Y, features)
        res_bin=CA.extract_features_fdr(path+'_chi2_binary.txt',50000,direction=True,binarization=True, allow_subseq=False)
        res_med=CA.extract_features_fdr(path+'_chi2_median.txt',50000,direction=True,binarization='median', allow_subseq=False)
        pos_bin,neg_bin=self.extract_top_but_n(res_bin,50000)
        pos_med,neg_med=self.extract_top_but_n(res_med,50000)
        Seq2Tax.write_in_fastafile(path+'_chi2_binary.fasta',res_bin)
        Seq2Tax.write_in_fastafile(path+'_chi2_relative.fasta',res_med)
        Seq2Tax.write_in_file(path+'_chi2_binary_taxa.txt',pos_bin,neg_bin)
        Seq2Tax.write_in_file(path+'_chi2_median_taxa.txt',pos_med,neg_med)

    @staticmethod
    def write_in_fastafile(filename,res, min_length=50):
        corpus=[]
        labels=[]
        for seq, score, pval, _, _ in res:
            if len(seq)>min_length:
                corpus.append(seq)
                labels.append(' '.join(['+' if score>0 else '-','p-val:'+str(pval)]))
        FileUtility.create_fasta_file(filename,corpus,labels)

    @staticmethod
    def write_in_file(filename,pos,neg):
        lines=[['direction','marker','p-value','genus']]
        for marker, pval,genus in pos:
            lines.append(['+',marker, str(pval),genus])
        for marker, pval,genus in neg:
            lines.append(['-',marker, str(pval),genus])
        FileUtility.save_list(filename, ['\t'.join(line) for line in lines])

    def micro_seq_to_tax(self, seq):
        tax_dict=self.seq2tax[seq]
        result_dict=dict()
        for level,level_details in tax_dict.items():
                    pairs=[(pairs.split('###')[0],float(pairs.split('###')[1]))for pairs in level_details.split()]
                    for pair in pairs:
                        (tax, prob)=pair
                        if prob>0.8:
                            if level not in result_dict:
                                result_dict[level]=dict()
                            if tax not in result_dict[level]:
                                result_dict[level][tax]=0
                            result_dict[level][tax]+=prob
        if 'species' in result_dict:
            return result_dict['species']
        elif 'genus' in result_dict:
            return result_dict['genus']
        elif 'family' in result_dict:
            return result_dict['family']
        elif 'order' in result_dict:
            return result_dict['order']
        elif 'class' in result_dict:
            return result_dict['class']
        elif 'phylum' in result_dict:
            return result_dict['phylum']
        else:
            return -1

    def extract_top_but_n(self, chi2_res, n, threshold=3.8):
        # 5
        pos=[]
        neg=[]
        idx=0
        while len(pos)<n and idx < len(chi2_res):
            if chi2_res[idx][1]>threshold  and chi2_res[idx][0] in self.seq2tax:
                if 'genus' in self.seq2tax[chi2_res[idx][0]]:
                    dict_tax={pairs.split('###')[0]:float(pairs.split('###')[1]) for pairs in self.seq2tax[chi2_res[idx][0]]['genus'].split()}
                    if np.max(list(dict_tax.values())) > 0.95:
                        pos.append((chi2_res[idx][0], chi2_res[idx][2], get_max_of_dict(dict_tax)))
                if 'family' in self.seq2tax[chi2_res[idx][0]]:
                    dict_tax={pairs.split('###')[0]:float(pairs.split('###')[1]) for pairs in self.seq2tax[chi2_res[idx][0]]['family'].split()}
                    if np.max(list(dict_tax.values())) > 0.95:
                        pos.append((chi2_res[idx][0], chi2_res[idx][2], get_max_of_dict(dict_tax)))
            idx+=1
        idx=0
        while len(neg)<n and idx < len(chi2_res):
            if chi2_res[idx][1]<-threshold and idx< len(chi2_res) and chi2_res[idx][0] in self.seq2tax:
                if 'genus' in self.seq2tax[chi2_res[idx][0]]:
                    dict_tax={pairs.split('###')[0]:float(pairs.split('###')[1]) for pairs in self.seq2tax[chi2_res[idx][0]]['genus'].split()}
                    if np.max(list(dict_tax.values())) > 0.95:
                        neg.append((chi2_res[idx][0], chi2_res[idx][2], get_max_of_dict(dict_tax)))
                elif 'family' in self.seq2tax[chi2_res[idx][0]]:
                    dict_tax={pairs.split('###')[0]:float(pairs.split('###')[1]) for pairs in self.seq2tax[chi2_res[idx][0]]['family'].split()}
                    if np.max(list(dict_tax.values())) > 0.95:
                        neg.append((chi2_res[idx][0], chi2_res[idx][2], get_max_of_dict(dict_tax)))
            idx+=1
        return pos, neg


