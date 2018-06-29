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

    parser.add_argument('-r','--reverse', action='store_true', help="Invert selection\n")
    parser.add_argument('-c','--complement', action='store_true', help="Invert selection\n")

    return vars(parser.parse_args())

def revcomp( par ):

    with utils.openw( par['out_f'] ) as outf:
        if par['complement'] and par['reverse']:
            res = (r.reverse_complement(id=r.id,description="RC") for r in SeqIO.parse( utils.openr(par['inp_f']), "fasta"))
        elif par['reverse']:
            res = (r.reverse(id=r.id,description="R") for r in SeqIO.parse( utils.openr(par['inp_f']), "fasta"))
        elif par['complement']:
            res = (r.complement(id=r.id,description="C") for r in SeqIO.parse( utils.openr(par['inp_f']), "fasta"))
        else:
            res = []

        SeqIO.write(res, outf, "fasta")

if __name__ == '__main__':
    params = read_params(sys.argv)
    revcomp( params )

