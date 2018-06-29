#!/usr/bin/env python

import sys
import collections
import utils

try:
    import argparse as ap
    import bz2 
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Convert Usearch ".uc" files in tab-delimited'
            ' files with the seed as first field followed by the other IDs\n')

    p.add_argument( 'uc', nargs='?', default=None, type=str,
            help=   "the input uc file [stdin if not present]")
    p.add_argument('txt', nargs='?', default=None, type=str,
            help=   "the output txt file compressed if fiven with bz2 extension\n"
                    "[stdout if not present]")

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )

    #with (open(args['uc']) if args['uc'] else sys.stdin) as inp:
    with utils.openr(args['uc']) as inp:
        for type,cln,seql,pid,strand,ing1,ign2,aln,query,target in (l.split('\t') for l in inp):
            if type == 'H':
                uc2cl[target.strip()].add( query )
            elif type == 'S' and  query not in uc2cl:
                uc2cl[query] = set()

    #openw = bz2.BZ2File if args['txt'].endswith(".bz2") else open
    with utils.openw(args['txt']) as out:
        for k,v in sorted(uc2cl.items(),key=lambda x:-len(x[1])):
            out.write( "\t".join([k]+list(v)) +"\n" )
