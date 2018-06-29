#!/usr/bin/env python

import sys

try:
    import argparse as ap
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

try:
    import pyphlan as ppa
except ImportError:
    sys.stderr.write( "pyphlan.py not found\n" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Get the leaves of a tree')

    p.add_argument( 'intree', nargs='?', default=None, type=str,
            help=   "the input tree [stdin if not present]")
    p.add_argument('out_file', nargs='?', default=None, type=str,
            help=   "the out file [stdout if not present]")

    """
    p.add_argument('-f', default=None, type=str,
            help=   "file containing the list of taxonomic clades to test " 
                    "end the corresponding leaf taxa" )

    st = ['lca','ltcs']
    p.add_argument('-s', choices=st, default='lca', type=str,
            help=  "select the lcs strategy")
    """
    return vars( p.parse_args() )


if __name__ == "__main__":
    args = read_params( sys.argv )
    ppa = ppa.PpaTree( args['intree'] )
    res = ppa.get_subtree_leaves( args['intree'])
    with (open( args['out_file'], "w") if args['out_file'] else sys.stdout) as outf:
        for r in res:
            #outf.write( "\t".join([r[0]]+r[1]) +  "\n" ) 
            outf.write( "\t".join(r[1]) +  "\n" ) 

