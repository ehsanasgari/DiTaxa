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
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

def read_params(args):
    parser = argparse.ArgumentParser(description='')
    arg = parser.add_argument
    arg( 'inp_f', metavar='INPUT_FILE', nargs='?', default=None, type=str,
         help="the input fna file [stdin if not present]")
    arg( 'out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str,
         help="the output fna file [stdout if not present]")
    arg( '-a', default=None, type=int,
         help="number of char after the match to report")
    arg( '-n', default=None, type=int,
         help="number of matching primers")

    parser.add_argument('-s', metavar='Subsequene to look for', required = True, type = str )

    return vars(parser.parse_args())

if __name__ == '__main__':
    par = read_params(sys.argv)

    ss = par['s'].lower()
    ssr = Seq(par['s']).reverse_complement().lower()
    f = os.path.basename(par['inp_f']).split(".")[0]
    with utils.openw( par['out_f'] ) as outf:
        for r in SeqIO.parse( utils.openr(par['inp_f']), "fasta"):
            rl = r.seq.lower()
            if ss in rl or ssr in rl:
                if par['a']:
                    if ss in rl:
                        i = str(rl).index(str(ss))
                        subs = rl[i:i+len(ss)+par['a']] if i+len(ss)+par['a'] < len(rl) else rl[i:]
                    else:
                        i = str(rl).index(str(ssr))
                        subs = rl[i:i+len(ssr)+par['a']] if i+len(ssr)+par['a'] < len(rl) else rl[i:]
                    outf.write( f + "\t" + str(r.id) + "\t" + str(subs) + "\n" )
                else:
                    if par['n']:
                        n = str(rl).count(str(ss)) + str(rl).count(str(ssr))
                        outf.write( f + "\t" + str(r.id) + "\t" + str(n) + "\n" )
                    else:
                        outf.write( f + "\t" + str(r.id) + "\n" )
