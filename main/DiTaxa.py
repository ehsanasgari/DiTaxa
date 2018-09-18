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
import warnings
import sys
sys.path.append('../')
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from utility.file_utility import FileUtility
#from utility.visualization_utility import plot_scatter
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from classifier.classical_classifiers import RFClassifier,SVM, LogRegression
#from classifier.DNN import DNN

class DiTaxaWorkflow:
    '''
        DiTaxaWorkflow
    '''

    def __init__(self, file_directory, file_extenstion, output_directory, dbname, vocab_size, seg_train_depth ,rep_sampling_depth, blastn_path, num_p=1,onlyfiles=[], override=1):
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
        print(str(len(self.fasta_files)), ' fasta files found in', self.file_directory)

        self.dbname=dbname
        self.vocab_size=vocab_size
        self.seg_train_depth=seg_train_depth
        self.rep_sampling_depth=rep_sampling_depth
        self.num_p=num_p
        self.output_directory=output_directory
        self.output_directory_inter=(output_directory[0:-1] if output_directory[-1]=='/' else output_directory)+'/intermediate_files/'
        self.blastn_path=blastn_path

        DiTaxaWorkflow.ensure_dir(self.output_directory)
        if not os.path.exists(self.output_directory+'logfile.txt'):
            self.log_file=[]
        else:
            self.log_file=FileUtility.load_list(self.output_directory+'logfile.txt')
        print('\t✔ DiTaxa workflow is getting started')

    def train_npe(self):
        '''
        :return:
        '''
        if self.override==1 or not DiTaxaWorkflow.exists(self.output_directory_inter+'npe_segmentatation/'):
            print('\t✔ Segmentation inference started.. ')
            start = time.time()
            G16s = NPESegmentTrainMetagenomics(self.file_directory, self.file_extenstion)
            DiTaxaWorkflow.ensure_dir(self.output_directory_inter+'npe_segmentatation/')
            G16s.generate(self.vocab_size, self.seg_train_depth,
                          self.output_directory_inter+'npe_segmentatation/'+self.dbname+'_'+'_'.join(['unique',str(self.vocab_size),'v',str(self.seg_train_depth),'s']),
                          backend='Sentencepiece',num_p=self.num_p)
            end = time.time()
            spent = (end - start)
            self.log_file.append('Segmentation inference '+'_'.join(['unique',str(self.vocab_size),'v',str(self.seg_train_depth),'s '])+str(spent)+' seconds , using '+str(self.num_p)+' cores')
        else:
            print('\t✔ Segmentation results directory exists. Thus, the step was bypassed')
            self.log_file.append('Segmentation results directory exists. Thus, the step was bypassed')
        FileUtility.save_list(self.output_directory+'logfile.txt',self.log_file)

    def representation_npe(self):
        '''
        :return:
        '''
        if self.override==1 or not DiTaxaWorkflow.exists(self.output_directory_inter+'npe_representation/'):
            print('\t✔ Creating NPE representations ...')
            start = time.time()
            G16s = NPESegmentApplyMetagenomics(self.file_directory, self.file_extenstion,self.output_directory_inter+'npe_segmentatation/'+self.dbname+'_'+'_'.join(['unique',str(self.vocab_size),'v',str(self.seg_train_depth),'s.model']),sampling_number=self.rep_sampling_depth,num_p=self.num_p)
            DiTaxaWorkflow.ensure_dir(self.output_directory_inter+'npe_representation/')
            G16s.generate_npes_all(save=self.output_directory_inter+'npe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth))
            end = time.time()
            spent = end-start
            print('\t✔ Generating the NPE representations at npe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth)+'  '+str(spent)+' seconds , using '+str(self.num_p)+'cores')
            self.log_file.append('Generating the NPE representations at npe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth)+'  '+str(spent)+' seconds , using '+str(self.num_p)+'cores')
        else:
            print('\t✔ Representation are already created. Thus, this is step is skipped!')
            self.log_file.append('Representation are already created. Thus, this is step is skipped!')
        FileUtility.save_list(self.output_directory+'logfile.txt',self.log_file)
        DiTaxaWorkflow.temp_cleanup()

    def biomarker_extraction(self, labeler, label_mapper, phenoname, p_value_threshold=0.05, pos_label=None, neg_label=None, excel=0):
        '''

        :return:
        '''
        print('\t✔ NPE Marker detection is started..')
        start = time.time()
        rep_base_path=self.output_directory_inter+'npe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth)
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
        DiTaxaWorkflow.ensure_dir(self.output_directory_inter+'npe_marker_files/')

        if self.override==1 or not DiTaxaWorkflow.exists(self.output_directory_inter+'npe_marker_files/'+'_'.join([phenoname,'chi2_relative.fasta'])):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                G16s = NPEMarkerDetection(rep_base_path +'.npz', rep_base_path +'_' + phenoname + '_Y.txt', rep_base_path + '_features', self.output_directory_inter + 'npe_marker_files/' + phenoname, selected_samples)
                G16s.extract_markers()

            end = time.time()
            spent = end-start
            print('\t✔ biomarker extraction ' + phenoname + '  ' + str(spent) + ' seconds , using ' + str(self.num_p) + ' cores')
            self.log_file.append('biomarker extraction ' + phenoname + '  ' + str(spent) + ' seconds , using ' + str(self.num_p) + ' cores')
        else:
            print('\t✔ Biomarker are already extracted. Thus, the statistical test was bypassed')
            self.log_file.append(' Biomarker are already extracted. Thus, the statistical test was bypassed')

        FileUtility.save_list(self.output_directory+'logfile.txt',self.log_file)

        print('\t✔ Taxonomic assignment of the markers..')

        if callable(labeler):
            phenotypes=[labeler(filenames[sample_id]) for sample_id in selected_samples]
        else:
            phenotypes=[labeler[filenames[sample_id]] for sample_id in selected_samples]

        fasta_file=self.output_directory_inter+'npe_marker_files/' + phenoname + '_chi2_relative.fasta'
        matrix_path=rep_base_path+'.npz'
        feature_file_path=rep_base_path+'_features'

        if len(FileUtility.read_fasta_sequences(fasta_file))>2000:
            remove_redundants=False
        else:
            remove_redundants=True

        FileUtility.ensure_dir(self.output_directory+'final_outputs/save_states/')
        if self.override==1 or not DiTaxaWorkflow.exists(self.output_directory+'final_outputs/save_states/'+phenoname+'.pickle'):
            start = time.time()
            Final_OBJ=NPEMarkerAnlaysis(fasta_file, matrix_path, feature_file_path, phenotypes, label_mapper, selected_samples, p_value_threshold=p_value_threshold, remove_redundants=remove_redundants,num_p=self.num_p, blastn_path=self.blastn_path)
            end = time.time()
            spent = end-start
            DiTaxaWorkflow.ensure_dir(self.output_directory+'final_outputs/')
            FileUtility.save_obj(self.output_directory+'final_outputs/save_states/' + phenoname, Final_OBJ)
            print('\t✔ Marker analysis and alignment ' + phenoname + '  ' + str(spent) + ' seconds, using ' + str(self.num_p) + 'cores')
            self.log_file.append('Marker analysis and alignment ' + phenoname + '  ' + str(spent) + ' seconds, using ' + str(self.num_p) + 'cores')
        else:
            Final_OBJ=FileUtility.load_obj(self.output_directory+'final_outputs/save_states/'+phenoname+'.pickle')
            print('\t✔ The aligned markers already existed and are loaded!')
            self.log_file.append('The aligned markers already existed and are loaded!')
        FileUtility.save_list(self.output_directory+'logfile.txt',self.log_file)

        # generating the tree
        Final_OBJ.generate_tree(self.output_directory +'final_outputs/', phenoname)

        if excel==1:
            print('\t✔ Creating marker excel file..')
            Final_OBJ.generate_excel(self.output_directory +'final_outputs/' + phenoname + '.xlsx',phenoname)
            X_addr=self.output_directory_inter+'npe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth)+'.npz'
            feature_addr=self.output_directory_inter+'npe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth)+'_features'
            markers=self.output_directory_inter+'npe_marker_files/'+phenoname+'_finalmarker_list.txt'
            Y=self.output_directory_inter+'npe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth)+'_'+phenoname+"_Y.txt"
            print ('\t✔ Creating t-sne plot..')
            DiTaxaWorkflow.plot_res(self.output_directory +'final_outputs/' + phenoname + '_tsne.pdf',X_addr,feature_addr,markers,Y, labels=['Negative', 'Positive'])

        if pos_label and neg_label:
            print('\t✔ Creating marker heatmap..')
            Final_OBJ.update_matrix_by_markers_N()
            Final_OBJ.generate_heatmap(self.output_directory +'final_outputs/' + phenoname + '_heatmap', pos_label=pos_label, neg_label=neg_label)
            if not excel ==1:
                print ('\t✔ Creating t-sne plot..')
                DiTaxaWorkflow.plot_res(self.output_directory +'final_outputs/' + phenoname + '_tsne.pdf',X_addr,feature_addr,markers,Y, labels=[neg_label, pos_label])
        DiTaxaWorkflow.temp_cleanup()
        print('\t⬛ Marker detection and analysis completed. You can find the results at '+self.output_directory+', in partuclar at final_outputs subdirectory.')

    @staticmethod
    def plot_scatter(ax, X, Y, x_label, y_label, title,legend_hide=True, legend_loc=4, label_dict=False, legend_size=7, legend_col=1, color_schemes_idx=1):
        global color_schemes
        color_schemes=[['green','blue','red','gold', 'cyan'], ['#ff0505', '#f2a041', '#cdff05', '#04d9cb', '#45a8ff', '#8503a6', '#590202', '#734d02', '#4ab304', '#025359', '#0454cc', '#ff45da', '#993829', '#ffda45', '#1c661c', '#05cdff', '#1c2f66', '#731f57', '#b24a04', '#778003', '#0e3322', '#024566', '#0404d9', '#e5057d', '#66391c', '#31330e', '#3ee697', '#2d7da6', '#20024d', '#33011c']+list(({'aliceblue':            '#F0F8FF','antiquewhite':         '#FAEBD7','aqua':                 '#00FFFF','aquamarine':           '#7FFFD4','azure':                '#F0FFFF','beige':                '#F5F5DC','bisque':               '#FFE4C4','black':                '#000000','blanchedalmond':       '#FFEBCD','blue':                 '#0000FF','blueviolet':           '#8A2BE2','brown':                '#A52A2A','burlywood':            '#DEB887','cadetblue':            '#5F9EA0','chartreuse':           '#7FFF00','chocolate':            '#D2691E','coral':                '#FF7F50','cornflowerblue':       '#6495ED','cornsilk':             '#FFF8DC','crimson':              '#DC143C','cyan':                 '#00FFFF','darkblue':             '#00008B','darkcyan':             '#008B8B','darkgoldenrod':        '#B8860B','darkgray':             '#A9A9A9','darkgreen':            '#006400','darkkhaki':            '#BDB76B','darkmagenta':          '#8B008B','darkolivegreen':       '#556B2F','darkorange':           '#FF8C00','darkorchid':           '#9932CC','darkred':              '#8B0000','darksalmon':           '#E9967A','darkseagreen':         '#8FBC8F','darkslateblue':        '#483D8B','darkslategray':        '#2F4F4F','darkturquoise':        '#00CED1','darkviolet':           '#9400D3','deeppink':             '#FF1493','deepskyblue':          '#00BFFF','dimgray':              '#696969','dodgerblue':           '#1E90FF','firebrick':            '#B22222','floralwhite':          '#FFFAF0','forestgreen':          '#228B22','fuchsia':              '#FF00FF','gainsboro':            '#DCDCDC','ghostwhite':           '#F8F8FF','gold':                 '#FFD700','goldenrod':            '#DAA520','gray':                 '#808080','green':                '#008000','greenyellow':          '#ADFF2F','honeydew':             '#F0FFF0','hotpink':              '#FF69B4','indianred':            '#CD5C5C','indigo':               '#4B0082','ivory':                '#FFFFF0','khaki':                '#F0E68C','lavender':             '#E6E6FA','lavenderblush':        '#FFF0F5','lawngreen':            '#7CFC00','lemonchiffon':         '#FFFACD','lightblue':            '#ADD8E6','lightcoral':           '#F08080','lightcyan':            '#E0FFFF','lightgoldenrodyellow': '#FAFAD2','lightgreen':           '#90EE90','lightgray':            '#D3D3D3','lightpink':            '#FFB6C1','lightsalmon':          '#FFA07A','lightseagreen':        '#20B2AA','lightskyblue':         '#87CEFA','lightslategray':       '#778899','lightsteelblue':       '#B0C4DE','lightyellow':          '#FFFFE0','lime':                 '#00FF00','limegreen':            '#32CD32','linen':                '#FAF0E6','magenta':              '#FF00FF','maroon':               '#800000','mediumaquamarine':     '#66CDAA','mediumblue':           '#0000CD','mediumorchid':         '#BA55D3','mediumpurple':         '#9370DB','mediumseagreen':       '#3CB371','mediumslateblue':      '#7B68EE','mediumspringgreen':    '#00FA9A','mediumturquoise':      '#48D1CC','mediumvioletred':      '#C71585','midnightblue':         '#191970','mintcream':            '#F5FFFA','mistyrose':            '#FFE4E1','moccasin':             '#FFE4B5','navajowhite':          '#FFDEAD','navy':                 '#000080','oldlace':              '#FDF5E6','olive':                '#808000','olivedrab':            '#6B8E23','orange':               '#FFA500','orangered':            '#FF4500','orchid':               '#DA70D6','palegoldenrod':        '#EEE8AA','palegreen':            '#98FB98','paleturquoise':        '#AFEEEE','palevioletred':        '#DB7093','papayawhip':           '#FFEFD5','peachpuff':            '#FFDAB9','peru':                 '#CD853F','pink':                 '#FFC0CB','plum':                 '#DDA0DD','powderblue':           '#B0E0E6','purple':               '#800080','red':                  '#FF0000','rosybrown':            '#BC8F8F','royalblue':            '#4169E1','saddlebrown':          '#8B4513','salmon':               '#FA8072','sandybrown':           '#FAA460','seagreen':             '#2E8B57','seashell':             '#FFF5EE','sienna':               '#A0522D','silver':               '#C0C0C0','skyblue':              '#87CEEB','slateblue':            '#6A5ACD','slategray':            '#708090','snow':                 '#FFFAFA','springgreen':          '#00FF7F','steelblue':            '#4682B4','tan':                  '#D2B48C','teal':                 '#008080','thistle':              '#D8BFD8','tomato':               '#FF6347','turquoise':            '#40E0D0','violet':               '#EE82EE','wheat':                '#F5DEB3','white':                '#FFFFFF','whitesmoke':           '#F5F5F5','yellow':               '#FFFF00','yellowgreen':          '#9ACD32'}).keys()),['#ff0505', '#f2a041', '#cdff05', '#04d9cb', '#45a8ff', '#8503a6', '#590202', '#734d02', '#4ab304', '#025359', '#0454cc', '#ff45da', '#993829', '#ffda45', '#1c661c', '#05cdff', '#1c2f66', '#731f57', '#b24a04', '#778003', '#0e3322', '#024566', '#0404d9', '#e5057d', '#66391c', '#31330e', '#3ee697', '#2d7da6', '#20024d', '#33011c']]
        #plt.rc('text', usetex=True)
        matplotlib.rcParams['mathtext.fontset'] = 'stix'
        matplotlib.rcParams['font.family'] = 'STIXGeneral'
        matplotlib.rcParams['mathtext.fontset'] = 'custom'
        matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
        matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
        matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
        matplotlib.rcParams["axes.edgecolor"] = "black"
        matplotlib.rcParams["axes.linewidth"] = 0.6
        matplotlib.rcParams["axes.titlesize"] = 20
        matplotlib.rcParams["axes.labelsize"] = 20
        matplotlib.rcParams["figure.facecolor"] = 'white'
        target=list(set(Y))
        target.sort()
        color_idx=[target.index(x) for x in Y]
        color_list=color_schemes[color_schemes_idx]

        for current_color in range(len(target)):
            color=color_list
            current_idxs=[idx for idx,v in enumerate(color_idx) if v==current_color]
            if label_dict:
                ax.scatter(X[current_idxs, 0], X[current_idxs, 1], c=color[current_color], label=label_dict[target[current_color]], cmap='viridis', alpha=0.4, edgecolors=None)
            else:
                ax.scatter(X[current_idxs, 0], X[current_idxs, 1], c=color[current_color], label=target[current_color], cmap='viridis', alpha=0.4, edgecolors=None)
        ax.set_title(title)
        rect = ax.patch
        rect.set_facecolor('white')
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_xticks([])
        ax.set_yticks([])


    @staticmethod
    def get_tsne(X, random_seed=0):
        X_tsne = TSNE(n_components=2, perplexity=40, verbose=0, learning_rate=10,random_state=random_seed).fit_transform(X.toarray())
        return X_tsne

    @staticmethod
    def plot_res(file_address, X_addr,features_addr, selected_addr, label_addr, labels=['Negative','Positive']):
        global color_schemes
        color_schemes=[['green','blue','red','gold', 'cyan'], ['#ff0505', '#f2a041', '#cdff05', '#04d9cb', '#45a8ff', '#8503a6', '#590202', '#734d02', '#4ab304', '#025359', '#0454cc', '#ff45da', '#993829', '#ffda45', '#1c661c', '#05cdff', '#1c2f66', '#731f57', '#b24a04', '#778003', '#0e3322', '#024566', '#0404d9', '#e5057d', '#66391c', '#31330e', '#3ee697', '#2d7da6', '#20024d', '#33011c']+list(({'aliceblue':            '#F0F8FF','antiquewhite':         '#FAEBD7','aqua':                 '#00FFFF','aquamarine':           '#7FFFD4','azure':                '#F0FFFF','beige':                '#F5F5DC','bisque':               '#FFE4C4','black':                '#000000','blanchedalmond':       '#FFEBCD','blue':                 '#0000FF','blueviolet':           '#8A2BE2','brown':                '#A52A2A','burlywood':            '#DEB887','cadetblue':            '#5F9EA0','chartreuse':           '#7FFF00','chocolate':            '#D2691E','coral':                '#FF7F50','cornflowerblue':       '#6495ED','cornsilk':             '#FFF8DC','crimson':              '#DC143C','cyan':                 '#00FFFF','darkblue':             '#00008B','darkcyan':             '#008B8B','darkgoldenrod':        '#B8860B','darkgray':             '#A9A9A9','darkgreen':            '#006400','darkkhaki':            '#BDB76B','darkmagenta':          '#8B008B','darkolivegreen':       '#556B2F','darkorange':           '#FF8C00','darkorchid':           '#9932CC','darkred':              '#8B0000','darksalmon':           '#E9967A','darkseagreen':         '#8FBC8F','darkslateblue':        '#483D8B','darkslategray':        '#2F4F4F','darkturquoise':        '#00CED1','darkviolet':           '#9400D3','deeppink':             '#FF1493','deepskyblue':          '#00BFFF','dimgray':              '#696969','dodgerblue':           '#1E90FF','firebrick':            '#B22222','floralwhite':          '#FFFAF0','forestgreen':          '#228B22','fuchsia':              '#FF00FF','gainsboro':            '#DCDCDC','ghostwhite':           '#F8F8FF','gold':                 '#FFD700','goldenrod':            '#DAA520','gray':                 '#808080','green':                '#008000','greenyellow':          '#ADFF2F','honeydew':             '#F0FFF0','hotpink':              '#FF69B4','indianred':            '#CD5C5C','indigo':               '#4B0082','ivory':                '#FFFFF0','khaki':                '#F0E68C','lavender':             '#E6E6FA','lavenderblush':        '#FFF0F5','lawngreen':            '#7CFC00','lemonchiffon':         '#FFFACD','lightblue':            '#ADD8E6','lightcoral':           '#F08080','lightcyan':            '#E0FFFF','lightgoldenrodyellow': '#FAFAD2','lightgreen':           '#90EE90','lightgray':            '#D3D3D3','lightpink':            '#FFB6C1','lightsalmon':          '#FFA07A','lightseagreen':        '#20B2AA','lightskyblue':         '#87CEFA','lightslategray':       '#778899','lightsteelblue':       '#B0C4DE','lightyellow':          '#FFFFE0','lime':                 '#00FF00','limegreen':            '#32CD32','linen':                '#FAF0E6','magenta':              '#FF00FF','maroon':               '#800000','mediumaquamarine':     '#66CDAA','mediumblue':           '#0000CD','mediumorchid':         '#BA55D3','mediumpurple':         '#9370DB','mediumseagreen':       '#3CB371','mediumslateblue':      '#7B68EE','mediumspringgreen':    '#00FA9A','mediumturquoise':      '#48D1CC','mediumvioletred':      '#C71585','midnightblue':         '#191970','mintcream':            '#F5FFFA','mistyrose':            '#FFE4E1','moccasin':             '#FFE4B5','navajowhite':          '#FFDEAD','navy':                 '#000080','oldlace':              '#FDF5E6','olive':                '#808000','olivedrab':            '#6B8E23','orange':               '#FFA500','orangered':            '#FF4500','orchid':               '#DA70D6','palegoldenrod':        '#EEE8AA','palegreen':            '#98FB98','paleturquoise':        '#AFEEEE','palevioletred':        '#DB7093','papayawhip':           '#FFEFD5','peachpuff':            '#FFDAB9','peru':                 '#CD853F','pink':                 '#FFC0CB','plum':                 '#DDA0DD','powderblue':           '#B0E0E6','purple':               '#800080','red':                  '#FF0000','rosybrown':            '#BC8F8F','royalblue':            '#4169E1','saddlebrown':          '#8B4513','salmon':               '#FA8072','sandybrown':           '#FAA460','seagreen':             '#2E8B57','seashell':             '#FFF5EE','sienna':               '#A0522D','silver':               '#C0C0C0','skyblue':              '#87CEEB','slateblue':            '#6A5ACD','slategray':            '#708090','snow':                 '#FFFAFA','springgreen':          '#00FF7F','steelblue':            '#4682B4','tan':                  '#D2B48C','teal':                 '#008080','thistle':              '#D8BFD8','tomato':               '#FF6347','turquoise':            '#40E0D0','violet':               '#EE82EE','wheat':                '#F5DEB3','white':                '#FFFFFF','whitesmoke':           '#F5F5F5','yellow':               '#FFFF00','yellowgreen':          '#9ACD32'}).keys()),['#ff0505', '#f2a041', '#cdff05', '#04d9cb', '#45a8ff', '#8503a6', '#590202', '#734d02', '#4ab304', '#025359', '#0454cc', '#ff45da', '#993829', '#ffda45', '#1c661c', '#05cdff', '#1c2f66', '#731f57', '#b24a04', '#778003', '#0e3322', '#024566', '#0404d9', '#e5057d', '#66391c', '#31330e', '#3ee697', '#2d7da6', '#20024d', '#33011c']]
        X=FileUtility.load_sparse_csr(X_addr)
        features=FileUtility.load_list(features_addr)
        features_selected=FileUtility.load_list(selected_addr)
        idx=[features.index(x) for  x in features_selected if x in features]
        X_selected=X[:,idx]
        Y=FileUtility.load_list(label_addr)
        X_tsne=DiTaxaWorkflow.get_tsne(X)
        X_red_tsne=DiTaxaWorkflow.get_tsne(X_selected)

        f = plt.figure(figsize=(16,8))
        ax1 = f.add_subplot(121)
        ax2 = f.add_subplot(122)
        DiTaxaWorkflow.plot_scatter(ax1, X_tsne, Y, 't-SNE 1', 't-SNE 0', '(i) t-SNE over NPE representations',legend_hide=False,  legend_loc=9,  legend_size=10,label_dict={'0':labels[0],'1':labels[1]}, color_schemes_idx=0)
        DiTaxaWorkflow.plot_scatter(ax2, X_red_tsne, Y, 't-SNE 1', 't-SNE 0', '(ii) t-SNE over selected markers',legend_hide=False,  legend_loc=9,  legend_size=10,label_dict={'0':labels[0],'1':labels[1]}, color_schemes_idx=0)
        plt.savefig(file_address)
        plt.close()

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
        sys.stdout = sys.__stdout__


    @staticmethod
    def enablePrint():
        sys.stdout = open(os.devnull, 'w')

    def classify_DNN(self, phenoname,model,gpu_id, batchsize, epochs):
        X=self.output_directory_inter+'npe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth)+'.npz'
        Y=self.output_directory_inter+'npe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth)+'_'+phenoname+'_Y.txt'
        DiTaxaWorkflow.ensure_dir(self.output_directory_inter+'classifications/NN/')
        out=self.output_directory_inter+'classifications/NN/'+phenoname
        DiTaxaWorkflow.DNN_classifier(out, X,Y,model,gpu_id,epochs,batchsize)

    def classify_classic(self, phenoname, model, cores):
        X=self.output_directory_inter+'npe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth)+'.npz'
        Y=self.output_directory_inter+'npe_representation/'+self.dbname+'_uniquepiece_'+str(self.rep_sampling_depth)+'_'+phenoname+'_Y.txt'
        DiTaxaWorkflow.ensure_dir(self.output_directory_inter+'classifications/'+model+'/')
        out=self.output_directory_inter+'classifications/'+model+'/'+phenoname
        DiTaxaWorkflow.classical_classifier(out, X,Y,model,cores)

    @staticmethod
    def classical_classifier(out_dir, X_file, Y_file, model, cores):
        #
        X=FileUtility.load_sparse_csr(X_file)
        # labels
        Y=[int (y) for y in FileUtility.load_list(Y_file)]

        if model=='RF':
            #### Random Forest classifier
            MRF = RFClassifier(X, Y)
            # results containing the best parameter, confusion matrix, best estimator, results on fold will be stored in this address
            MRF.tune_and_eval(out_dir, njobs=cores)
        elif model=='SVM':
            #### Support Vector Machine classifier
            MSVM = SVM(X, Y)
            # results containing the best parameter, confusion matrix, best estimator, results on fold will be stored in this address
            MSVM.tune_and_eval(out_dir, njobs=cores)
        elif model=='LR':
            #### Logistic regression classifier
            MLR = LogRegression(X, Y)
            # results containing the best parameter, confusion matrix, best estimator, results on fold will be stored in this address
            MLR.tune_and_eval(out_dir, njobs=cores)

    @staticmethod
    def DNN_classifier(out_dir, X_file, Y_file, arch, gpu_id, epochs, batch_size):
        # k-mer data
        X=FileUtility.load_sparse_csr(X_file).toarray()
        # labels
        Y=[int (y) for y in FileUtility.load_list(Y_file)]
        DeepNN=DNN(X,Y,model_arch=arch)
        DeepNN.cross_validation(out_dir, gpu_dev=gpu_id, n_fold=10, epochs=epochs, batch_size=batch_size, model_strct='mlp')

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
