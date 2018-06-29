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
    p = ap.ArgumentParser(description='Convert txt files to libsvm\n')

    p.add_argument( 'txt', nargs='?', default=None, type=str,
            help=   "the input txt file [stdin if not present]")
    p.add_argument('ls', nargs='?', default=None, type=str,
            help=   "the output ilibsvm file compressed if fiven with bz2 extension\n"
                    "[stdout if not present]")

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )

    with utils.openr(args['txt']) as inp:
        data = zip(*[l.strip().split('\t') for l in inp])
        outd = [[d[0]]+[str(i+1)+":"+dd for i,dd in enumerate(d[1:])] for d in data[1:]]
        with utils.openw(args['ls']) as out:
            for o in outd:
                out.write( "\t".join(o) +"\n" )
