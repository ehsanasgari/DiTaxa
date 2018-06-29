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
    p = ap.ArgumentParser(description='Convert core gene txt file'
            ' substituting gene IDs with genomes IDs\n')

    p.add_argument( 't2g', nargs='?', default=None, type=str,
            help=   "")
    p.add_argument('g2t', nargs='?', default=None, type=str,
            help=   "")

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )

    g2t = {}

    with utils.openr( args['t2g'] ) as inp:
        for ll in (l.strip().split('\t') for l in inp):
            to = int(ll[0])
            for g in ll[1:]:
                g2t[int(g)] = to 
    
    with utils.openw(args['g2t']) as out:
        for g,t in g2t.iteritems():
            out.write( "\t".join([str(g),str(t)]) +"\n" )
