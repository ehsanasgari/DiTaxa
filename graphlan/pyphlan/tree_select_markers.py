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
            description='Select top markers')

    p.add_argument( 'infile', nargs='?', default=None, type=str,
            help=   "the input marker file [stdin if not present]")
    p.add_argument('outfile', nargs='?', default=None, type=str,
            help=   "the output core file [stdout if not present]")
    p.add_argument('-n', metavar="Maximum number of markers to be selected"
                                 "for each clade",
            default=100, type=int )
    p.add_argument('--th', metavar="Threshold on markerness",
            default=None, type=str )

    return vars( p.parse_args() )


if __name__ == "__main__":
    args = read_params( sys.argv )
    tree = ppa.PpaTree( None )
    markers = tree.select_markers( args['infile'], markerness_th = args['th'], max_markers = args['n'] )

    with open( args['outfile'], "w" ) as outf:
        for clade in markers:
            for m in clade:
                outf.write( "\n".join( m ) + "\n" )


