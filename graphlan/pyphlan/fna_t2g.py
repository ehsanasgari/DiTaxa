#!/usr/bin/env python 
from __future__ import with_statement

import sys
import argparse
import os
import textwrap
from collections import namedtuple as nt
import random as rnd
rnd.seed(1982)
import utils
from Bio import SeqIO

def read_params(args):
    parser = argparse.ArgumentParser(description='List the genes in the genome file')
    arg = parser.add_argument
    arg( 'inp_f', metavar='INPUT_FILE', default=None, type=str,
         help="the input fna file")
    arg( 'out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str,
         help="the output txt file [stdout if not present]")

    return vars(parser.parse_args())

def genome_id( fn ):
    return str(-int(os.path.basename(fn).split(".")[0]))

if __name__ == '__main__':
    par = read_params(sys.argv)
  
    
    ids = [r.id for r in SeqIO.parse( utils.openr(par['inp_f']), "fasta")]

    with utils.openw( par['out_f']) as out:
        out.write( "\t".join( [genome_id(par['inp_f'])]+ids ) + "\n" )

