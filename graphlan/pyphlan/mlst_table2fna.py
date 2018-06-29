#!/usr/bin/env python

import sys
import collections
import utils
try:
    import argparse as ap
    import random
    from Bio import SeqIO 
    from Bio.SeqRecord import SeqRecord
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Create a fasta file with the'
            'concatenated mlst sequence from a mlst table and the single sequences')

    p.add_argument( '--fna', required=True, default=None, type=str,
            help=   "the file with all the MLST profiles [in the format >profilineName_profileID")
    p.add_argument( '--txt', required=True, default=None, type=str,
            help=   "the table of the samples to profiles [tab-delimited, columns ID are profileName]")
    p.add_argument( '--nmiss', default = 0, type = int )

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )

    fna = SeqIO.to_dict(SeqIO.parse( utils.openr(args['fna']), "fasta"))
    fna_out = [] 

    profiles = {}
    mlst_names = []
    with utils.openr(args['txt']) as inp:
        for i,line in enumerate(inp):
            if i == 0:
                mlst_names = line.strip().split('\t')[1:]
                continue
            l = line.strip().split('\t')
            profiles[l[0]] = dict([(na,l[n+1]) for n,na in enumerate(mlst_names)])
    for s,p in profiles.items():
        seq = ""
        skip = 0 
        for n in mlst_names:
            name = p[n]
            if name not in fna:
                name = n+"_"+p[n]
            if name in fna:
                seq += fna[name].seq
            else:
                skip += 1
                continue
        if skip > args['nmiss']:
            continue
        sample = SeqRecord( "sample_"+s  )
        sample.id = "t"+s if s[0] in "0123456789" else s
        sample.description = ""
        sample.seq = seq
        fna_out.append( sample )
    
    SeqIO.write(fna_out, sys.stdout, "fasta")
