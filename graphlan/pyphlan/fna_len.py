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
from Bio.SeqRecord import SeqRecord

def read_params(args):
    parser = argparse.ArgumentParser(description='Display the length of each fasta entry')
    arg = parser.add_argument
    arg( 'inp_f', metavar='INPUT_FILE', nargs='?', default=None, type=str,
         help="the input fna file [stdin if not present]")
    arg( 'out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str,
         help="the output txt file [stdout if not present]")
    arg( '-t','--total', action='store_true', help="Print only the sum of the length of all sequences\n") 
    return vars(parser.parse_args())

if __name__ == '__main__':
    par = read_params(sys.argv)

    ltot = 0
    with utils.openw( par['out_f'] ) as outf:
        for r in SeqIO.parse( utils.openr(par['inp_f']), "fasta"):
            l = len(r.seq)
            if par['total']:
                ltot += l
            else:
                outf.write( "\t".join( [r.id,str(l)] ) + "\n" ) 
        if par['total']:
            outf.write( str(ltot) + "\n" ) 
