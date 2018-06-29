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
    sys.stderr.write( "pyphlan.py not found" )
    sys.exit(-1)


def read_params( args ):
    p = ap.ArgumentParser(
            description='Screen core genes for markers')

    p.add_argument( 'intree', nargs='?', default=None, type=str,
            help=   "the input tree [stdin if not present]")
    p.add_argument('outfile', nargs='?', default=None, type=str,
            help=   "the output core file [stdout if not present]")
    p.add_argument('--hitmap', metavar="File containing global hits",
            default=None, type=str )
    p.add_argument('--cores', metavar="File containing unique cores",
            default=None, type=str )
    p.add_argument('--core_info', metavar="File containing unique core info",
            default=None, type=str )

    return vars( p.parse_args() )


if __name__ == "__main__":
    args = read_params( sys.argv )
    tree = ppa.PpaTree( args['intree'] )
    cores = tree.find_markers(args['cores'], args['hitmap'], args['core_info'])

    with open( args['outfile'], "w" ) as outf:
        for k,v in sorted(cores.items(),key=lambda x:x[0]):
            if len(v) > 1:
                continue
            outf.write( "\t".join( [str(k)] + [str(s) for s in v]) + "\n" )
