__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu"
__project__ = "LLP - DiTaxa"
__website__ = "https://llp.berkeley.edu/ditaxa/"

import argparse
import os
import os.path
import sys
from main.DiTaxa import DiTaxaWorkflow
from utility.file_utility import FileUtility


def checkArgs(args):
    '''
        This function checks the input arguments and returns the errors (if exist) otherwise reads the parameters
    '''
    # keep all errors
    err = "";
    # Using the argument parser in case of -h or wrong usage the correct argument usage
    # will be prompted
    parser = argparse.ArgumentParser()

    ## to do : chi2 print

    # input directory #################################################################################################
    parser.add_argument('--indir', action='store', dest='input_dir', default=False, type=str,
                        help='directory of 16S rRNA samples')

    # file type #######################################################################################################
    parser.add_argument('--ext', action='store', dest='filetype', default='fastq', type=str,
                        help='extension of the sample files, the default is fastq')

    # to override the previous files or to continue ####################################################################
    parser.add_argument('--override', action='store', dest='override',default=1, type=int,
                        help='Override the existing files?')

    # output directory #################################################################################################
    parser.add_argument('--outdir', action='store', dest='output_dir', default=False, type=str,
                        help="directory for storing the output files, if doesn't exist will be created.")

    # dbname ################################################################################################
    parser.add_argument('--dbname', action='store', dest='dbname', default=False, type=str,
                        help='dataset name: to be used for figures and output creation!')

    # cores ################################################################################################
    parser.add_argument('--cores', action='store', dest='cores', default=4, type=int,
                        help='Number of cores to be used, default is 4')

    # label filename #################################################################################################
    parser.add_argument('--fast2label', action='store', dest='fast2label', default=False, type=str,
                        help='tabular mapping between fatsa/fastq file names and their labels')

    # blast path #################################################################################################
    parser.add_argument('--blastn', action='store', dest='blastn', default=False, type=str,
                        help='path to the bin directory of blastn; get the latest from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/')

    # label values ##################################################################################
    parser.add_argument('--phenomap', action='store', dest='phenomap', default=None, type=str,
                        help='pair of comma separated label:[0 or 1]. e.g., untreated_disease:1, treated_diesease:0, healthy:0, ..')

    # label values ##################################################################################
    parser.add_argument('--phenoname', action='store', dest='phenoname', default=None, type=str,
                        help='Phenotype setting name, if not given the labeling scheme will be used.')

    # generate heatmap ###############################################################################
    parser.add_argument('--heatmap', action='store', dest='heatmap', default=None, type=str,
                        help='to generate heatmap of top 100 markers: positive_title:negative_title')

    # generate excel file  ###############################################################################
    parser.add_argument('--excel', action='store', dest='excel', default=1, type=int,
                        help='to generate excel output')


    ######################################################################################################
    parser.add_argument('--classify', action='store', dest='classify', type=str, default=False,
                    choices=[False, 'RF', 'SVM', 'DNN', 'LR'],
                    help='train_predictor: choice of classifier from RF, SVM, DNN')

    parser.add_argument('--batchsize', action='store', dest='batch_size', type=int, default=10,
                        help='train_predictor-model/DNN: batch size for deep learning')

    parser.add_argument('--gpu_id', action='store', dest='gpu_id', type=str, default='0',
                        help='train_predictor-model/DNN: GPU id for deep learning')

    parser.add_argument('--epochs', action='store', dest='epochs', type=int, default=100,
                        help='train_predictor-model/DNN: number of epochs for deep learning')

    parser.add_argument('--arch', action='store', dest='dnn_arch', type=str, default='1024,0.2,512',
                        help='train_predictor-model/DNN: The comma separated definition of neural network layers connected to eahc other, you do not need to specify the input and output layers, values between 0 and 1 will be considered as dropouts')


    parsedArgs = parser.parse_args()

    if (not os.access(parsedArgs.input_dir, os.F_OK)):
        err = err + "\nError: Permission denied or could not find the directory!"
        return err
    if (not os.access(parsedArgs.fast2label, os.F_OK)):
        err = err + "\nError: Permission to the label file is denied!"
        return err
    try:
        label_dict = dict()
        pheno_temp=[]
        for x in parsedArgs.phenomap.split(','):
            k, n = x.split(':')
            k = k
            n = int(n)
            label_dict[k]=n
            pheno_temp.append('@'.join([k,str(n)]))
        pheno_temp='#'.join(pheno_temp)
        if not parsedArgs.phenoname:
            phenoname=pheno_temp
        else:
            phenoname=parsedArgs.phenoname
    except:
        err = err + "\nWrong format for labels!"
        return err
    if parsedArgs.heatmap:
        if not len(parsedArgs.heatmap.split(':'))==2:
            return err + "\nThe heatmap inputs is incorrect!"

    if not os.access(parsedArgs.blastn, os.F_OK):
        print('The blast path is incorrect..')
        exit()

    Pipeline = DiTaxaWorkflow(parsedArgs.input_dir,
                                          parsedArgs.filetype,parsedArgs.output_dir,parsedArgs.dbname, 50000,5000,-1, parsedArgs.blastn, num_p=parsedArgs.cores, override=parsedArgs.override)
    Pipeline.train_npe()
    Pipeline.representation_npe()
    labels={line.split()[0].split('/')[-1]:line.split()[1] for line in FileUtility.load_list(parsedArgs.fast2label)}

    if parsedArgs.heatmap:
        pos_label, neg_label =parsedArgs.heatmap.split(':')
        Pipeline.biomarker_extraction(labels,label_dict,phenoname, excel=parsedArgs.excel, pos_label=pos_label,neg_label=neg_label)
    else:
        print (parsedArgs.excel)
        exit()
        Pipeline.biomarker_extraction(labels,label_dict,phenoname, excel=parsedArgs.excel)

    if parsedArgs.classify:
        print('Classification requested..')
        if parsedArgs.classify=='DNN':
            '''
                Deep learning
            '''
            arch=[int(layer) if float(layer)>1 else float(layer) for layer in parsedArgs.dnn_arch.split(',')]
            Pipeline.classify_DNN(phenoname, arch, parsedArgs.gpu_id,parsedArgs.batch_size,parsedArgs.epochs)
        else:
            '''
                SVM and Random Forest
            '''
            if parsedArgs.classify in ['SVM','RF','LR']:
                Pipeline.classify_classic(phenoname, parsedArgs.classify, parsedArgs.cores)
            else:
                return  "\nNot able to recognize the model!"

if __name__ == '__main__':
    err = checkArgs(sys.argv)
    if err:
        print(err)
        exit()
