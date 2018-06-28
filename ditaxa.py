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


def checkArgs(args):
    '''
        This function checks the input arguments and returns the errors (if exist) otherwise reads the parameters
    '''
    # keep all errors
    err = "";
    # Using the argument parser in case of -h or wrong usage the correct argument usage
    # will be prompted
    parser = argparse.ArgumentParser()

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
    parser.add_argument('--fastlist', action='store', dest='filelist', default=False, type=str,
                        help='fasta file names, the same order should be used in the label file')

    # label file #################################################################################################
    parser.add_argument('--label', action='store', dest='label', default=False, type=str,
                        help='labels with the same order as the files in --fastlist')

    # label values ##################################################################################
    parser.add_argument('--label_vals', action='store', dest='label_value', default=None, type=str,
                        help='pair of comma separated label:[0 or 1]. e.g., untreated_disease:1, treated_diesease:0, healthy:0, ..')


    parsedArgs = parser.parse_args()

    if (not os.access(parsedArgs.input_dir, os.F_OK)):
        err = err + "\nError: Permission denied or could not find the directory!"
        return err
    if (not os.access(parsedArgs.label, os.F_OK)):
        err = err + "\nError: Permission to the label file is denied!"
        return err
    if (not os.access(parsedArgs.filelist, os.F_OK)):
        err = err + "\nError: Permission to the label file is denied!"
        return err
    try:
        label_dict = dict()
        for x in parsedArgs.label_value.split(','):
            k, n = x.split(':')
            k = k
            n = int(n)
            label_dict[k]=n
    except:
        err = err + "\nWrong format for labels!"
        return err

    Pipeline = DiTaxaWorkflow(parsedArgs.input_dir,
                                          parsedArgs.filetype,parsedArgs.output_dir,parsedArgs.dbname,50000,5000,-1,num_p=parsedArgs.cores, override=parsedArgs.override)
    Pipeline.train_npe()
    Pipeline.representation_npe()
    labels=FileUtility.load_list(parsedArgs.label)
    labels={x.split('/')[-1]:labels[idx] for idx,x in enumerate(FileUtility.load_list(parsedArgs.filelist))}
    Pipeline.biomarker_extraction(labels,label_dict,parsedArgs.dbname)


if __name__ == '__main__':
    err = checkArgs(sys.argv)
    if err:
        print(err)
        exit()

