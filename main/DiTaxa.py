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
from make_representations.npe_segmentation_train import NPESegmentTrainMetagenomics
from make_representations.npe_segmentation_apply import NPESegmentApplyMetagenomics
from marker_detection.npe_biomarker_extraction import NPEMarkerDetection
from marker_detection.npe_generate_taxa_tree import NPEMarkerAnlaysis
import sys, os
import shutil


class DiTaxaWorkflow:
    '''
        DiTaxaWorkflow
    '''

    def __init__(self, file_directory, file_extenstion, output_directory, dbname, vocab_size, seg_train_depth ,rep_sampling_depth, num_p=1,onlyfiles=[], override=1):
        '''
        :param file_directory: the samples directory
        :param file_extenstion: the file extension fastq or fasta
        :param onlyfiles: filter a list of files
        :param backend: which backend to use
        '''
        self.override=override
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

        DiTaxaWorkflow.ensure_dir(self.output_directory)
        if not os.path.exists(self.output_directory+'logfile.txt'):
            self.log_file=[]
        else:
            self.log_file=FileUtility.load_list(self.output_directory+'logfile.txt')
        print('DiTaxa workflow is getting started')

    def train_npe(self):
        '''
        :return:
        '''
        if self.override==1 or not DiTaxaWorkflow.exists(self.output_directory+'npe_segmentatation/'):
            print('npe training started.. ')
            DiTaxaWorkflow.blockPrint()
            start = time.time()
            G16s = NPESegmentTrainMetagenomics(self.file_directory, self.file_extenstion)
            DiTaxaWorkflow.ensure_dir(self.output_directory+'npe_segmentatation/')
            G16s.generate(self.vocab_size, self.seg_train_depth,
                          self.output_directory+'npe_segmentatation/'+self.dbname+'_'+'_'.join(['unique',str(self.vocab_size),'v',str(self.seg_train_depth),'s']),
                          backend='Sentencepiece',num_p=self.num_p)
            end = time.time()
            spent = (end - start)
            self.log_file.append('training segmentation '+'_'.join(['unique',str(self.vocab_size),'v',str(self.seg_train_depth),'s '])+str(spent)+' seconds , using '+str(self.num_p)+' cores')
            DiTaxaWorkflow.enablePrint()
        else:
            print('segmentation directory already exists and the training was bypassed')
            self.log_file.append('segmentation directory already exists and the training was bypassed')
        FileUtility.save_list(self.output_directory+'logfile.txt',self.log_file)

    def representation_npe(self):
        '''
        :return:
        '''
        if self.override==1 or not DiTaxaWorkflow.exists(self.output_directory+'npe_representation/'):
            print('npe generation started..')
            start = time.time()
            G16s = NPESegmentApplyMetagenomics(self.file_directory, self.file_extenstion,self.output_directory+'npe_segmentatation/'+self.dbname+'_'+'_'.join(['unique',str(self.vocab_size),'v',str(self.seg_train_depth),'s.model']),sampling_number=self.rep_sampling_depth,num_p=self.num_p)
            DiTaxaWorkflow.ensure_dir(self.output_directory+'npe_representation/')
            G16s.generate_npes_all(save=self.output_directory+'npe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth))
            end = time.time()
            spent = end-start
            self.log_file.append('generating the representations npe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth)+'  '+str(spent)+' seconds , using '+str(self.num_p)+'cores')
        else:
            print('representation directory already exists and this was bypassed')
            self.log_file.append('representation directory already exists and the step was bypassed')
        FileUtility.save_list(self.output_directory+'logfile.txt',self.log_file)
        DiTaxaWorkflow.temp_cleanup()

    def biomarker_extraction(self, labeler, label_mapper, phenoname, p_value_threshold=0.05, pos_label=None, neg_label=None):
        '''

        :return:
        '''
        print('npe marker detection started')
        start = time.time()
        rep_base_path=self.output_directory+'npe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth)
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

        FileUtility.save_list(rep_base_path +'_' + phenoname + '_Y.txt', Y)
        DiTaxaWorkflow.ensure_dir(self.output_directory+'npe_marker_files/')

        if self.override==1 or not DiTaxaWorkflow.exists(self.output_directory+'npe_marker_files/'+'_'.join([phenoname,'chi2_relative.fasta'])):
            DiTaxaWorkflow.blockPrint()
            G16s = NPEMarkerDetection(rep_base_path +'.npz', rep_base_path +'_' + phenoname + '_Y.txt', rep_base_path + '_features', self.output_directory + 'npe_marker_files/' + phenoname, selected_samples)
            G16s.extract_markers()
            DiTaxaWorkflow.enablePrint()
            end = time.time()
            spent = end-start
            self.log_file.append('biomarker extraction ' + phenoname + '  ' + str(spent) + ' seconds , using ' + str(self.num_p) + ' cores')
        else:
            print('Biomarker file already exists and the statistical test was bypassed')
            self.log_file.append('Biomarker file already exists and the statistical test was bypassed')

        FileUtility.save_list(self.output_directory+'logfile.txt',self.log_file)

        print('npe marker taxonomic detection is getting started..')

        if callable(labeler):
            phenotypes=[labeler(filenames[sample_id]) for sample_id in selected_samples]
        else:
            phenotypes=[labeler[filenames[sample_id]] for sample_id in selected_samples]

        fasta_file=self.output_directory+'npe_marker_files/' + phenoname + '_chi2_relative.fasta'
        matrix_path=rep_base_path+'.npz'
        feature_file_path=rep_base_path+'_features'

        if len(FileUtility.read_fasta_sequences(fasta_file))>2000:
            remove_redundants=False
        else:
            remove_redundants=True

        if self.override==1 or not DiTaxaWorkflow.exists(self.output_directory+'final_outputs/'+phenoname+'.pickle'):
            start = time.time()
            Final_OBJ=NPEMarkerAnlaysis(fasta_file, matrix_path, feature_file_path, phenotypes, label_mapper, selected_samples, p_value_threshold=p_value_threshold, remove_redundants=remove_redundants,num_p=self.num_p)
            end = time.time()
            spent = end-start
            DiTaxaWorkflow.ensure_dir(self.output_directory+'final_outputs/')
            FileUtility.save_obj(self.output_directory +'final_outputs/' + phenoname, Final_OBJ)
            self.log_file.append('blasting extraction ' + phenoname + '  ' + str(spent) + ' seconds, using ' + str(self.num_p) + 'cores')
        else:
            Final_OBJ=FileUtility.load_obj(self.output_directory+'final_outputs/'+phenoname+'.pickle')
            print('The aligned markers already existed and are loaded!')
            self.log_file.append('The aligned markers already existed and are loaded!')
        FileUtility.save_list(self.output_directory+'logfile.txt',self.log_file)

        # generating the tree
        Final_OBJ.generate_tree(self.output_directory +'final_outputs/', phenoname)



        if pos_label and neg_label:
            print('Creating marker heatmap..')
            Final_OBJ.update_matrix_by_markers_N()
            Final_OBJ.generate_heatmap(self.output_directory +'final_outputs/' + phenoname + '_heatmap', pos_label=pos_label, neg_label=neg_label)
        DiTaxaWorkflow.temp_cleanup()


    @staticmethod
    def temp_cleanup():
        for the_file in os.listdir('tmp/'):
            file_path = os.path.join('tmp/', the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(e)

    @staticmethod
    def exists(file_path):
        return os.path.exists(file_path)

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
    Pipeline = DiTaxaWorkflow('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/dental/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/dentaloutput_new/','periodontal',50000,5000,-1,num_p=20)
    Pipeline.train_npe()
    Pipeline.representation_npe()
    f=(lambda x: 'Periodental' if 'd' in x else 'Healthy')
    Pipeline.biomarker_extraction(f,{'Periodental':1,'Healthy':0},'Periodontal', p_value_threshold=0.05, pos_label='Periodontal',neg_label='Healthy')

def dental_tree():
    Pipeline = DiTaxaWorkflow('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/dental/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/dentaloutput/','periodontal',50000,5000,-1,num_p=20)
    #Pipeline.train_npe()
    #Pipeline.representation_npe()
    f=(lambda x: 'Periodental' if 'd' in x else 'Healthy')
    Pipeline.biomarker_extraction(f,{'Periodental':1,'Healthy':0},'Periodontal')

def RA():
    Pipeline = DiTaxaWorkflow('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/RA/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/RAoutput/','RA',50000,5000,-1,num_p=20)
    Pipeline.train_npe()
    Pipeline.representation_npe()
    labels=FileUtility.load_list('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/ra/rep/labels.txt')
    labels={x.split('/')[-1]:labels[idx] for idx,x in enumerate(FileUtility.load_list('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/ra/rep/ra_selfposnpe_10000_npe_5000_meta'))}
    Pipeline.biomarker_extraction(labels,{'untreated_RA':1,'healthy':0,'treated_RA':0,'psoriatic':0},'RA_vs_all')


def RA_healthy():
    Pipeline = DiTaxaWorkflow('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/RA/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/RAoutput/','RA',50000,5000,-1,num_p=20)
    #Pipeline.train_npe()
    #Pipeline.representation_npe()
    labels=FileUtility.load_list('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/ra/rep/labels.txt')
    labels={x.split('/')[-1]:labels[idx] for idx,x in enumerate(FileUtility.load_list('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/ra/rep/ra_selfposnpe_10000_npe_5000_meta'))}
    Pipeline.biomarker_extraction(labels,{'untreated_RA':1,'treated_RA':0},'untreated_vs_treated')


def IBD():
    Pipeline = DiTaxaWorkflow('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/crohn/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/IBDout/','IBD',50000,5000,-1,num_p=20)
    Pipeline.train_npe()
    Pipeline.representation_npe()
    labels=dict([(x.split()[0]+'.fastq',x.split()[1]) for x in FileUtility.load_list('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/crohns/rep/Crohns_lables.txt')])
    Pipeline.biomarker_extraction(labels,{'CD':1,'no':0,'control':0},'CD_vs_healthy')

def IBD_rep():
    Pipeline = DiTaxaWorkflow('/mounts/data/proj/asgari/dissertation/datasets/deepbio/microbiome/crohn/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/IBDout/','IBD',50000,5000,-1,num_p=20)
    #Pipeline.train_npe()
    #Pipeline.representation_npe()
    labels=dict([(x.split()[0]+'.fastq',x.split()[1]) for x in FileUtility.load_list('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/crohns/rep/Crohns_lables.txt')])
    Pipeline.biomarker_extraction(labels,{'IC':1,'no':0,'control':0},'IC_vs_healthy')

def synthetic():
    Pipeline = DiTaxaWorkflow('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/synthetic/16S_evalutaion_1000/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/SyntheticOUT/','synthetic',50000,5000,-1,num_p=20)
    Pipeline.train_npe()
    Pipeline.representation_npe()
    f=(lambda x: 'Case' if 'case' in x else 'control')
    Pipeline.biomarker_extraction(f,{'case':1,'control':0},'case_vs_control')

def synthetic_test():
    Pipeline = DiTaxaWorkflow('/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/synthetic/16S_evalutaion_1000/',
                                          'fastq','/mounts/data/proj/asgari/dissertation/git_repos/16S_datasets/SyntheticOUT/','synthetic',50000,5000,-1,num_p=20)
    #Pipeline.train_npe()
    #Pipeline.representation_npe()
    f=(lambda x: 'case' if 'case' in x else 'control')
    Pipeline.biomarker_extraction(f,{'case':1,'control':0},'case_vs_control')

if __name__ == '__main__':
    dental()
