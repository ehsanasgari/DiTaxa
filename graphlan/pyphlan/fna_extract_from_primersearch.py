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


def read_params( args ):
    parser = argparse.ArgumentParser(description='Amplify DNA regions using primersearch output\n')

    arg = parser.add_argument 
    arg( '--ps', required = True, default=None, type=str )
    arg( 'fna', nargs='?', default=None, type=str,
            help=   "the input fna file [stdin if not present]")
    arg('out', nargs='?', default=None, type=str,
            help=   "the output fna file\n"
                    "[stdout if not present]")

    return vars( parser.parse_args() )


def parse_primersearch( fn ):
    seqs = {}
    cur,seq,fs,rs,al,rseq,fseq = None, None, None, None, None, None, None
    with utils.openr( fn, "U" ) as inpf:
        for l in inpf:
            line = l.strip()
            if line.startswith("Amplimer") and 'Amplimer length' not in line:
                cur = line
            elif line.startswith("Sequence"):
                seq = line.split("Sequence:")[1].strip()
            elif 'hits forward strand at ' in line:
                fseq = line.split()[0]
                fs = int(line.split("hits forward strand at ")[1].split("with")[0])
            elif 'hits reverse strand at ' in line:
                rseq = line.split()[0]
                rs = int(line.split("hits reverse strand at ")[1].split("with")[0].strip()[1:-1])
            elif 'Amplimer length' in line:
                al = int(line.split("Amplimer length: ")[1].split("bp")[0])
                seqs[cur] = { 'seq' : seq, 'fs' : fs, 'rs' : rs, 'al' : al, 'rseq' : rseq, 'fseq' : fseq }
                cur,seq,fs,rs,al,rseq,fseq = None, None, None, None, None, None, None
    return seqs


if __name__ == "__main__":
    args = read_params( sys.argv )

    extr = parse_primersearch( args['ps']  )
    
    seqs2extr = {}
    for k,v in extr.items():
        if v['seq'] in seqs2extr:
            seqs2extr[v['seq']][k] = v
        else:
            seqs2extr[v['seq']] = { k: v }

    with utils.openw( args['out'] ) as outf:
        for r in SeqIO.parse( utils.openr(args['fna']), "fasta"):
            if r.id in seqs2extr:
                for pn,ext in seqs2extr[r.id].items():
                    sq = SeqRecord( r.id )
                    sq.id = r.id + " " + pn
                    sq.description = r.description + " " + pn
                    sq.seq = r.seq[ ext['fs']+len(ext['fseq']):len(r.seq)-ext['rs']-len(ext['rseq'])]
                    SeqIO.write(sq, outf, "fasta") 





