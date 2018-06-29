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
import regex

def read_params(args):
    parser = argparse.ArgumentParser(description='')
    arg = parser.add_argument
    arg( 'inp_f', metavar='INPUT_FILE', nargs='?', default=None, type=str,
         help="the input fna file [stdin if not present]")
    arg( 'out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str,
         help="the output fna file [stdout if not present]")
    parser.add_argument('-l', metavar='', required = True, type = str )
    parser.add_argument('-r', metavar='', required = True, type = str )
    parser.add_argument('-e', metavar='', required = False, type = int, default = 0 )

    return vars(parser.parse_args())

if __name__ == '__main__':
    par = read_params(sys.argv)

    rn,ln = "", ""
    if par['r'].count(":"):
        rn,par['r'] = par['r'].split(":")
        rn+=":"
    if par['l'].count(":"):
        ln,par['l'] = par['l'].split(":")
        ln+=":"

    ne = str(par['e'])
    c_r = par['r'].lower()
    c_r_rev = Seq(par['r']).reverse_complement().lower()
    c_l = par['l'].lower()
    c_l_rev = Seq(par['l']).reverse_complement().lower()
    f = os.path.basename(par['inp_f']).split(".")[0]

    with utils.openw( par['out_f'] ) as outf:
        for seq in SeqIO.parse( utils.openr(par['inp_f']), "fasta"):
            seql = str(seq.seq.lower())
            r_rev,l_rev,r,l = c_r_rev,c_l_rev,c_r,c_l      
           
            rr = regex.findall( "("+r+"){e<="+ne+"}", seql )
            lr = regex.findall( "("+l+"){e<="+ne+"}", seql )
            r_revr = regex.findall( "("+str(r_rev)+"){e<="+ne+"}", seql )
            l_revr = regex.findall( "("+str(l_rev)+"){e<="+ne+"}", seql )

            if len(rr) > 1:
                outf.write( str(rr) +" unspecific 1\n" )
            if len(lr) > 1:
                outf.write( str(lr) +" unspecific 2\n" )
            if len(r_revr) > 1:
                outf.write( str(r_revr) +" unspecific 3\n" )
            if len(l_revr) > 1:
                outf.write( str(l_revr) +" unspecific 4\n" )
            if len(rr) == 1 and len(l_revr) == 1:
                lr,r_revr = l_revr, rr
            if len(lr) == 1 and len(r_revr) == 1:
                i = str(seql).index(str(lr[0])) + len(l)
                j = str(seql).index(str(r_revr[0])) - len(r_rev)
                if i > j:
                    i,j = j,i
                outf.write( ln+str(l) + "\t" + rn+str(r_rev) + "\t" + str(len(seql[i:j])) + "\t" + str(len(seql[i:j])+len(l)+len(r_rev)) +"\t" + str(seql[i:j]) +   "\n" )


    
    """
    with utils.openw( par['out_f'] ) as outf:
        for seq in SeqIO.parse( utils.openr(par['inp_f']), "fasta"):
            r_rev,l_rev,r,l = c_r_rev,c_l_rev,c_r,c_l      
            while True:
                seql = seq.seq.lower()
                if r in seql and l_rev in seql:
                    r_rev,l = r,l_rev
                if r_rev in seql and l in seql:
                    if seql.count( r_rev ) > 1 or seql.count( l ) > 1:
                        outf.write( "ambiguous primers\n" )
                        break
                    i = str(seql).index(str(l)) + len(l)
                    j = str(seql).index(str(r_rev)) - len(r_rev)
                    outf.write( ln+str(l) + "\t" + rn+str(r_rev) + "\t" + str(len(seql[i:j])) + "\t" + str(len(seql[i:j])+len(l)+len(r_rev)) +"\t" + str(seql[i:j]) +   "\n" )
                    break
                r_rev,l_rev = r_rev[:-1],l_rev[1:]
                l,r = l[1:],l[:1]
    """
