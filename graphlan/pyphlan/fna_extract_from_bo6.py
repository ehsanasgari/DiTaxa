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
    parser = argparse.ArgumentParser(description='Split/Select/Randomize/Subsample a multi fasta file')
    arg = parser.add_argument
    arg( 'inp_f', metavar='INPUT_FILE', nargs='?', default=None, type=str,
         help="the input fna file [stdin if not present]")
    arg( 'out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str,
         help="the output fna file [stdout if not present]")

    parser.add_argument('--extract_targets', action='store_true', help="Select fna entries\n")
    parser.add_argument('-i', action='store_true', help="Add hit stats to fna entries\n")
    
    parser.add_argument('--bo6', metavar='Bo6 file', required=True, type = str )

    return vars(parser.parse_args())

if __name__ == '__main__':
    par = read_params(sys.argv)

    inp_mat = (l.rstrip('\n').split("\t") for l in (utils.openr(par['bo6'])))

    if par['extract_targets']:
        toextr = ((l[1], l[2], l[3], l[11], int(l[8]), int(l[9])) for l in inp_mat)
    else:
        toextr = ((l[0], l[2],  l[3], l[11], int(l[6]), int(l[7])) for l in inp_mat)
   
    inpfasta = SeqIO.to_dict(SeqIO.parse( utils.openr(par['inp_f']), "fasta"))

    out_seqs = []
    for n,pid,l,bit,fr,to in toextr:
        n = inpfasta[n][min(fr,to):max(fr,to)]
        if par['i']:
            p = "_pid"+pid.strip()+"_l"+l.strip()+"_bs"+bit.strip()
        else:
            p = ""
        n.id = n.id+"_"+str(fr)+"_"+str(to)+p
        out_seqs.append( n )

    SeqIO.write(out_seqs, utils.openw(par['out_f']), "fasta") 



