#!/usr/bin/env python

import sys
import utils
import bz2 

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
    p = ap.ArgumentParser(description='Finds pairwise distances between nodes in the tree based on branch lengths')

    p.add_argument( 'intree', nargs='?', default=None, type=str,
            help=   "the input tree")
    p.add_argument('out_file', nargs='?', default=None, type=str,
            help=   "the output file (b2zipped if ending with '.bz2')\n"
                    "[stdout if not present]")
    p.add_argument( '-n',  action='store_true', 
                    help = "Distances normalized with respect to the total branch length" )

    return vars( p.parse_args() )


if __name__ == "__main__":
    args = read_params( sys.argv )
    ppatree = ppa.PpaTree( args['intree'] )

    dists = ppa.dist_matrix(ppatree.tree) 
    tbl = ppatree.tree.total_branch_length() if args['n'] else 1.0
    #tbl = ppatree.tree.total_branch_length()-1.0 if args['n'] else 1.0

    with utils.openw( args['out_file'] ) as out:
        for k1,v1 in dists.items():
            for k2,v2 in v1.items():
                if k1 < k2:
                    out.write( "\t".join([k1,k2,str(v2/tbl)]) +"\n" )    

