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

    Pipeline = DiTaxaWorkflow(parsedArgs.input_dir,
                                          parsedArgs.filetype,parsedArgs.output_dir,parsedArgs.dbname,50000,5000,-1,num_p=parsedArgs.cores, override=parsedArgs.override)
    Pipeline.train_npe()
    Pipeline.representation_npe()
    labels={line.split()[0].split('/')[-1]:line.split()[1] for line in FileUtility.load_list(parsedArgs.fast2label)}

    if parsedArgs.heatmap:
        pos_label, neg_label =parsedArgs.heatmap.split(':')
        Pipeline.biomarker_extraction(labels,label_dict,phenoname, excel=parsedArgs.excel, pos_label=pos_label,neg_label=neg_label)
    else:
        Pipeline.biomarker_extraction(labels,label_dict,phenoname, excel=parsedArgs.excel)


if __name__ == '__main__':
    err = checkArgs(sys.argv)
    if err:
        print(err)
        exit()

