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
    parser = argparse.ArgumentParser(description='Split/Select/Randomize/Subsample a multi fasta file')
    arg = parser.add_argument
    arg( 'inp_f', metavar='INPUT_FILE', nargs='?', default=None, type=str,
         help="the input fna file [stdin if not present]")
    arg( 'out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str,
         help="the output fna file [stdout if not present]")

    parser.add_argument('-l', metavar='Length of the prodiced reads', default=100, type = int )
    parser.add_argument('-s', metavar='Step in sampling reads', default=20, type = int )

    return vars(parser.parse_args())

if __name__ == '__main__':
    par = read_params(sys.argv)

    step = par['s']
    rlen = par['l']

    with utils.openw( par['out_f'] ) as outf:
        for r in SeqIO.parse( utils.openr(par['inp_f']), "fasta"):
            gene = r.seq
            oid = r.id
            reads = []
            for i in xrange(0, len(gene), step ):
                if len(gene) - i < rlen - step:
                    break
                nid = "_".join([oid,str(i)])
                reads.append( SeqRecord( gene[i:i+rlen], id = nid, description = "") )

            SeqIO.write(reads, outf, "fasta")
