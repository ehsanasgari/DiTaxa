__author__ = "Ehsaneddin Asgari"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu or ehsaneddin.asgari@helmholtz-hzi.de"
__project__ = "LLP - Life Language Processing"
__website__ = "https://llp.berkeley.edu/"

import sys
sys.path.append('../')
from utility.list_set_util import get_max_of_dict
import numpy as np
from chi2analysis.chi2analysis import Chi2Analysis
import tqdm
from multiprocessing import Pool
from Bio.Blast.Applications import NcbiblastnCommandline
from utility.file_utility import FileUtility
from Bio.Blast import NCBIXML
from scipy.sparse import csr_matrix
import os


class FastaTaxa(object):

    def __init__(self):
        self.ez_taxa_dict={x.split()[0]:x.split()[1].split(';') for x in FileUtility.load_list('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/EZ/raw/eztaxon_id_taxonomy.txt')}



    def extract_markers(self, path, X, Y, features,dataset, setting_des, min_length=50, top_n_features=50000,top=5):
        '''
        :param path:
        :param X:
        :param Y:
        :param features:
        :return:
        '''
        CA=Chi2Analysis(X, Y, features)
        res_bin=CA.extract_features_fdr(path+'_chi2_binary.txt',top_n_features,direction=True,binarization=True, allow_subseq=False)
        res_med=CA.extract_features_fdr(path+'_chi2_median.txt',top_n_features,direction=True,binarization='median', allow_subseq=False)
        pos_table, neg_table=self.table_top_but_n(res_bin,top)
        pos_med_table, neg_med_table=self.table_top_but_n(res_med,top)

        self.make_latex_table(path+'_binary_table' ,dataset, setting_des, pos_table, neg_table)
        self.make_latex_table(path+'_median_table', dataset, setting_des, pos_med_table, neg_med_table)
        FastaTaxa.write_in_fastafile(path+'_chi2_binary.fasta',res_bin,min_length)
        FastaTaxa.write_in_fastafile(path+'_chi2_relative.fasta',res_med,min_length)


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

    def table_top_but_n(self, chi2_res, n, max_level=100):
        final_table_pos=[]
        final_table_neg=[]
        i_pos=0
        i_neg=0
        for idx, (seq, xi2score, pval, _, _) in tqdm.tqdm(enumerate(chi2_res),total=max_level):
            if (xi2score >0 and i_pos<n) or (xi2score <0 and i_neg<n):
                FileUtility.create_fasta_file('temp.fasta',[seq],['temp'])
                blastx_cline=NcbiblastnCommandline(query='temp.fasta', db="/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/EZ/raw/eztaxon_qiime_full.fasta", evalue=0.001, outfmt=5, out="temp.xml")
                blastx_cline()
                f=open("temp.xml",'r')
                blast_records = NCBIXML.parse(f)
                flag=False
                score=-1
                alignment_length=-1
                results=[]
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            if not flag and score==-1:
                                score=hsp.score
                                alignment_length=hsp.align_length
                                flag=True

                            if hsp.score >= score and hsp.align_length>=alignment_length:
                                results.append((self.ez_taxa_dict[alignment.hit_id],hsp.expect))
                if len(results)>0:
                    res=FastaTaxa.is_not_too_ambiguous(results)
                    if res:
                        if xi2score>0:
                            if i_pos<n:
                                final_table_pos.append([seq,'+', pval,res])
                            i_pos+=1
                        else:
                            if i_neg<n:
                                final_table_neg.append([seq,'-', pval,res])
                            i_neg+=1
                    if (i_pos >=n or i_neg>=n) and idx > max_level:
                        return final_table_pos,final_table_neg
        return final_table_pos,final_table_neg

    @staticmethod
    def is_not_too_ambiguous(results):
        levels_id={'Superkingdom':1,'phylum':1,'class':2,'order':3,'family':4,'genus':5,'species':6}

        if len(set([x[0][levels_id['family']] for x in results]))>1:
            return False
        genuses=set([x[0][levels_id['genus']] for x in results])
        if len(genuses)==1:
            species=set([x[0][levels_id['species']] for x in results])
            if len(species)==1:
                return ';'.join(results[0][0])
            elif len(species)<10:
                return '\makecell[l]{' +';'.join(results[0][0][0:6])+'\\\\ \\textbf{Species: }\\small{'+', '.join(list(set([x[0][6] for x in results])))+'}}'
            else:
                return '\makecell[l]{' +';'.join(results[0][0][0:6])+'}'
        else:
            return False
        #elif len(genuses)<3:
        #    return '\makecell{' +';'.join(results[0][0][0:5])+'\\\\'+'\\\\'.join(list(set([';'.join(x[0][5:7]) for x in results])))+'}'
        #else:
        #    return ';'.join(results[0][0][0:5])

    def make_latex_table(self, path, dataset, setting_des, pos, neg):

            f=open(path,'w')
            f.write(('\multicolumn{1}{|c|}{\multirow{'+str(len(pos)+len(neg))+'}{*}{\makecell{'+dataset+' \\\\'+setting_des+'}}} &'+'\n').replace('_','$\-$'))
            for idx,(seq,sign, pval,res) in enumerate(pos):
                seq='\makecell[l]{'+('\\\\'.join([seq[i:i+60] for i in range(0, len(seq), 60)]))+'}'
                pval=round(pval,5)
                f.write((('&' if idx>0 else '')+' & '.join([str(sign), seq, str(pval), res])+'  \\\\ \cline{2-5}\n').replace('_','$\-$'))

            for idx, (seq,sign, pval,res) in enumerate(neg):
                pval=round(pval,5)
                seq='\makecell[l]{'+('\\\\'.join([seq[i:i+60] for i in range(0, len(seq), 60)]))+'}'
                f.write((' & '+' & '.join([str(sign), seq, str(pval), res])+'  \\\\  \cline{2-5}\n').replace('_','$\-$'))
            f.write('\\hline')
            f.close()

def get_matrix_label(X, labeler, SRX2row, pos=[], neg=[]):
    label_vec=[]
    X=X.toarray()
    New_X=[]
    for srx,label in labeler.items():
        if label in pos:
            New_X.append(X[SRX2row[srx],:])
            label_vec.append(1)
        elif label in neg:
            New_X.append(X[SRX2row[srx],:])
            label_vec.append(0)
    return csr_matrix(New_X), label_vec

import pandas as pd

def get_matrix_label_cr(X, labels, pos=[], neg=[]):
    label_vec=[]
    X=X.toarray()
    New_X=[]
    for row in range(X.shape[0]):
        label=labels[row]
        if label in pos:
            New_X.append(X[row,:])
            label_vec.append(1)
        elif label in neg:
            New_X.append(X[row,:])
            label_vec.append(0)
    return csr_matrix(New_X), label_vec

def RA():
    labeler=dict([tuple(x.split(',')) for x in FileUtility.load_list('../../bio_cpe_datacollection/representations/RA/grouping.csv')[1::]])
    SRX2row=dict([(x.split('/')[-1].split('.')[0],idx) for idx, x in enumerate(FileUtility.load_list('../../bio_cpe_datacollection/representations/RA/RA_100000_cpe_10000_meta'))])
    mat=FileUtility.load_sparse_csr('../../16S_datasets/ra/rep/ra_selfposcpe_50000_cpe_-1.npz')
    features=FileUtility.load_list('../../16S_datasets/ra/rep/ra_selfposcpe_50000_cpe_-1_features')
    X,labels=get_matrix_label(mat, labeler, SRX2row, pos=['treated_RA'], neg=['healthy'])
    SeqTaxa=FastaTaxa()
    SeqTaxa.extract_markers('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/ra/markers/tr_healthy',X,labels,features,'Rheumatoid Arthritis', 'treared vs. healthy')

def Crohns():
    mat=FileUtility.load_sparse_csr('../../16S_datasets/crohns/rep/crohns_selfcpe_50000_cpe_-1.npz')
    features=FileUtility.load_list('../../16S_datasets/crohns/rep/crohns_selfcpe_50000_cpe_-1_features')
    df=pd.read_csv('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/crohn_disease/convert_value', delimiter='\t')
    df2=pd.read_csv('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/crohn_disease/mapping_files/30808_mapping_file.txt',delimiter='\t')
    mapper=dict([tuple(x.split()[::-1]) for x in FileUtility.load_list('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/crohn_disease/labels')[1::]])
    mapper_err2ers=pd.Series(df.secondary_sample_accession.values,index=df.run_accession).to_dict()
    mapper_ers2label=pd.Series(df2.diagnosis.values,index=df2.SampleID).to_dict()
    mapper_err2label=dict([(err,mapper_ers2label[mapper[ers]]) for err,ers in mapper_err2ers.items()])
    labels=[mapper_err2label[x.split('/')[-1].split('.')[0]] for x in FileUtility.load_list('../../16S_datasets/crohns/rep/crohns_selfcpe_50000_cpe_-1_meta')]
    X,Y=get_matrix_label_cr(mat, labels,pos=['IC'],neg=['no','control'])
    SeqTaxa=FastaTaxa()
    SeqTaxa.extract_markers('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/crohns/marker/IC_healthy',X,Y,features,"Indeterminate colitis (IC)", 'IC vs. healthy')

def dental():
    mat=FileUtility.load_sparse_csr('../../16S_datasets/dental/rep/dental_selfcpe_50000_cpe_-1.npz')
    features=FileUtility.load_list('../../16S_datasets/dental/rep/dent_selfcpe_50000_cpe_5000_features')
    labels=[1 if 'd' in x.split('/')[-1] else 0 for x in FileUtility.load_list('../../16S_datasets/dental/rep/dent_selfcpe_50000_cpe_5000_meta')]
    SeqTaxa=FastaTaxa()
    SeqTaxa.extract_markers('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/dental/markers/d_healthy',mat,labels,features,"Periodontal disease", 'Periodontal disease vs. healthy')

if __name__ == '__main__':
    import os
    path = '/mounts/data/proj/asgari/dissertation/datasets/deepbio/taxonomy/ncbi-blast-2.5.0+/bin/'
    os.environ['PATH'] += ':'+path
    #os.system("export PATH=$PATH:/mounts/data/proj/asgari/dissertation/datasets/deepbio/taxonomy/ncbi-blast-2.5.0+/bin/")
    dental()
