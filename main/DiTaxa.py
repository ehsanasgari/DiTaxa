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
from make_representations.npe_efficient import train_npe
import sentencepiece as spm
import time
import os
from make_representations.biocpe_segmentation_train import NPESegmentTrainMetagenomics
from make_representations.biocpe_segmentation_apply import NPESegmentApplyMetagenomics
from marker_detection.biocpe_biomarker_extraction import NPEMarkerDetection
from marker_detection.biocpe_generate_taxa_tree import NPEMarkerAnlaysis
import sys, os



class DiTaxaWorkflow:
    '''
        DiTaxaWorkflow
    '''

    def __init__(self, file_directory, file_extenstion, output_directory, dbname, vocab_size, seg_train_depth ,rep_sampling_depth, num_p=1,onlyfiles=[]):
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

        self.dbname=dbname
        self.vocab_size=vocab_size
        self.seg_train_depth=seg_train_depth
        self.rep_sampling_depth=rep_sampling_depth
        self.num_p=num_p
        self.output_directory=output_directory

        Bio16SCPEPipeline.ensure_dir(self.output_directory)
        if not os.path.exists(self.output_directory+'logfile.txt'):
            self.log_file=[]
        else:
            self.log_file=FileUtility.load_list(self.output_directory+'logfile.txt')
        print('pipeline started')

    def train_cpe(self):
        '''
        :return:
        '''
        print('cpe training started.. it might take more than 1 hour for more than 1000 samples')
        Bio16SCPEPipeline.blockPrint()
        start = time.time()
        G16s = BioCPESegmentTrainMetagenomics(self.file_directory, self.file_extenstion)
        Bio16SCPEPipeline.ensure_dir(self.output_directory+'biocpe_segmentatation/')
        G16s.generate(self.vocab_size, self.seg_train_depth,
                      self.output_directory+'biocpe_segmentatation/'+self.dbname+'_'+'_'.join(['unique',str(self.vocab_size),'v',str(self.seg_train_depth),'s']),
                      backend='Sentencepiece',num_p=self.num_p)
        end = time.time()
        spent = (end - start)
        self.log_file.append('training segmentation '+'_'.join(['unique',str(self.vocab_size),'v',str(self.seg_train_depth),'s '])+str(spent)+' seconds , using '+str(self.num_p)+'cores')
        Bio16SCPEPipeline.enablePrint()
        FileUtility.save_list(self.output_directory+'logfile.txt',self.log_file)

    def representation_cpe(self):
        '''
        :return:
        '''
        print('cpe generation started..')
        start = time.time()
        G16s = BioCPESegmentApplyMetagenomics(self.file_directory, self.file_extenstion,self.output_directory+'biocpe_segmentatation/'+self.dbname+'_'+'_'.join(['unique',str(self.vocab_size),'v',str(self.seg_train_depth),'s.model']),sampling_number=self.rep_sampling_depth,num_p=self.num_p)
        Bio16SCPEPipeline.ensure_dir(self.output_directory+'biocpe_representation/')
        G16s.generate_cpes_all(save=self.output_directory+'biocpe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth))
        end = time.time()
        spent = end-start
        self.log_file.append('generating the representations biocpe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth)+'  '+str(spent)+' seconds , using '+str(self.num_p)+'cores')
        FileUtility.save_list(self.output_directory+'logfile.txt',self.log_file)


    def biomarker_extraction(self, labeler, label_mapper, name_setting,p_value_threshold=0.05, pos_label=None,neg_label=None):
        '''

        :return:
        '''
        print('cpe marker detection started')
        Bio16SCPEPipeline.blockPrint()
        start = time.time()
        rep_base_path=self.output_directory+'biocpe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth)
        filenames=[x.split('/')[-1] for x in FileUtility.load_list(rep_base_path+'_meta')]

        # CHECK EXISTING LABELS
        if callable(labeler):
            selected_samples=[idx for idx, file in enumerate(filenames) if labeler(file) in label_mapper]
        else:
            selected_samples=[idx for idx, file in enumerate(filenames) if labeler[file] in label_mapper]


        if callable(labeler):
            Y=[str(label_mapper[labeler(filenames[sample_id])]) for sample_id in selected_samples]
        else:
            Y=[str(label_mapper[labeler[filenames[sample_id]]]) for sample_id in selected_samples]

        FileUtility.save_list(rep_base_path+'_'+name_setting+'_Y.txt', Y)
        Bio16SCPEPipeline.ensure_dir(self.output_directory+'biocpe_marker_files/')
        G16s = BioCPEMarkerDetection(rep_base_path+'.npz',rep_base_path+'_'+name_setting+'_Y.txt',rep_base_path+'_features',self.output_directory+'biocpe_marker_files/'+name_setting, selected_samples)
        G16s.extract_markers()
        end = time.time()
        spent = end-start
        self.log_file.append('biomarker extraction '+name_setting+'  '+str(spent)+' seconds , using '+str(self.num_p)+'cores')
        FileUtility.save_list(self.output_directory+'logfile.txt',self.log_file)
        Bio16SCPEPipeline.enablePrint()
        print('cpe marker taxonomic detection started')
        start = time.time()

        if callable(labeler):
            phenotypes=[labeler(filenames[sample_id]) for sample_id in selected_samples]
        else:
            phenotypes=[labeler[filenames[sample_id]] for sample_id in selected_samples]

        fasta_file=self.output_directory+'biocpe_marker_files/'+name_setting+'_chi2_relative.fasta'
        matrix_path=rep_base_path+'.npz'
        feature_file_path=rep_base_path+'_features'

        if len(FileUtility.read_fasta_sequences(fasta_file))>2000:
            remove_redundants=False
        else:
            remove_redundants=True

        Final_OBJ=BioCPEMarkerAnlaysis(fasta_file, matrix_path, feature_file_path, phenotypes, label_mapper, selected_samples, p_value_threshold=p_value_threshold, remove_redundants=remove_redundants,num_p=self.num_p)
        end = time.time()
        spent = end-start
        Bio16SCPEPipeline.ensure_dir(self.output_directory+'final_outputs/')
        FileUtility.save_obj(self.output_directory+'final_outputs/'+name_setting,Final_OBJ)
        Final_OBJ.generate_tree(self.output_directory+'final_outputs/',name_setting)
        self.log_file.append('blasting extraction '+name_setting+'  '+str(spent)+' seconds, using '+str(self.num_p)+'cores')
        FileUtility.save_list(self.output_directory+'logfile.txt',self.log_file)
        if pos_label and neg_label:
            Final_OBJ.generate_heatmap(self.output_directory+'final_outputs/'+name_setting+'_heatmap', pos_label=pos_label, neg_label=neg_label)



    @staticmethod
    def ensure_dir(file_path):
        directory = os.path.dirname(file_path)
        if not os.path.exists(directory):
            os.makedirs(directory)
    @staticmethod
    def blockPrint():
        sys.stdout = open(os.devnull, 'w')

    @staticmethod
    def enablePrint():
        sys.stdout = sys.__stdout__


def dental():
    Pipeline = Bio16SCPEPipeline('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/dental/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/dentaloutput_new/','periodontal',50000,5000,-1,num_p=20)
    Pipeline.train_cpe()
    Pipeline.representation_cpe()
    f=(lambda x: 'Periodental' if 'd' in x else 'Healthy')
    Pipeline.biomarker_extraction(f,{'Periodental':1,'Healthy':0},'Periodontal', p_value_threshold=0.05, pos_label='Periodontal',neg_label='Healthy')

def dental_tree():
    Pipeline = Bio16SCPEPipeline('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/dental/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/dentaloutput/','periodontal',50000,5000,-1,num_p=20)
    #Pipeline.train_cpe()
    #Pipeline.representation_cpe()
    f=(lambda x: 'Periodental' if 'd' in x else 'Healthy')
    Pipeline.biomarker_extraction(f,{'Periodental':1,'Healthy':0},'Periodontal')

def RA():
    Pipeline = Bio16SCPEPipeline('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/RA/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/RAoutput/','RA',50000,5000,-1,num_p=20)
    Pipeline.train_cpe()
    Pipeline.representation_cpe()
    labels=FileUtility.load_list('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/ra/rep/labels.txt')
    labels={x.split('/')[-1]:labels[idx] for idx,x in enumerate(FileUtility.load_list('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/ra/rep/ra_selfposcpe_10000_cpe_5000_meta'))}
    Pipeline.biomarker_extraction(labels,{'untreated_RA':1,'healthy':0,'treated_RA':0,'psoriatic':0},'RA_vs_all')


def RA_healthy():
    Pipeline = Bio16SCPEPipeline('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/RA/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/RAoutput/','RA',50000,5000,-1,num_p=20)
    #Pipeline.train_cpe()
    #Pipeline.representation_cpe()
    labels=FileUtility.load_list('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/ra/rep/labels.txt')
    labels={x.split('/')[-1]:labels[idx] for idx,x in enumerate(FileUtility.load_list('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/ra/rep/ra_selfposcpe_10000_cpe_5000_meta'))}
    Pipeline.biomarker_extraction(labels,{'untreated_RA':1,'treated_RA':0},'untreated_vs_treated')


def IBD():
    Pipeline = Bio16SCPEPipeline('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/crohn/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/IBDout/','IBD',50000,5000,-1,num_p=20)
    Pipeline.train_cpe()
    Pipeline.representation_cpe()
    labels=dict([(x.split()[0]+'.fastq',x.split()[1]) for x in FileUtility.load_list('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/crohns/rep/Crohns_lables.txt')])
    Pipeline.biomarker_extraction(labels,{'CD':1,'no':0,'control':0},'CD_vs_healthy')

def IBD_rep():
    Pipeline = Bio16SCPEPipeline('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/crohn/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/IBDout/','IBD',50000,5000,-1,num_p=20)
    #Pipeline.train_cpe()
    #Pipeline.representation_cpe()
    labels=dict([(x.split()[0]+'.fastq',x.split()[1]) for x in FileUtility.load_list('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/crohns/rep/Crohns_lables.txt')])
    Pipeline.biomarker_extraction(labels,{'IC':1,'no':0,'control':0},'IC_vs_healthy')

def synthetic():
    Pipeline = Bio16SCPEPipeline('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/synthetic/16S_evalutaion_1000/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/SyntheticOUT/','synthetic',50000,5000,-1,num_p=20)
    Pipeline.train_cpe()
    Pipeline.representation_cpe()
    f=(lambda x: 'Case' if 'case' in x else 'control')
    Pipeline.biomarker_extraction(f,{'case':1,'control':0},'case_vs_control')

def synthetic_test():
    Pipeline = Bio16SCPEPipeline('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/synthetic/16S_evalutaion_1000/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/SyntheticOUT/','synthetic',50000,5000,-1,num_p=20)
    #Pipeline.train_cpe()
    #Pipeline.representation_cpe()
    f=(lambda x: 'case' if 'case' in x else 'control')
    Pipeline.biomarker_extraction(f,{'case':1,'control':0},'case_vs_control')

if __name__ == '__main__':
    dental()
