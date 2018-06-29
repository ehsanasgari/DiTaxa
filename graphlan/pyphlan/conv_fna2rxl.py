#!/usr/bin/env python

import sys
import collections
import utils
try:
    import argparse as ap
    import bz2 
    import random
    from Bio import SeqIO 
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Convert fasta files files in tab-delimited'
            ' files with the seed as first field followed by the other IDs\n')

    p.add_argument( 'fna', nargs='?', default=None, type=str,
            help=   "the input uc file [stdin if not present]")
    p.add_argument('rxl', nargs='?', default=None, type=str,
            help=   "the output txt file compresse if fiven with bz2 extension\n"
                    "[stdout if not present]")
    """
    p.add_argument('--subsample', metavar="Subsampling rate",
            default=1.0, type=float )
    p.add_argument('-n', metavar="Minimum number of matching taxa",
            default=0, type=int )
    p.add_argument('-p', metavar="Prefix for taxon names",
            default="", type=str )
    """
    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )


    fna = SeqIO.to_dict(SeqIO.parse( utils.openr(args['fna']), "fasta"))

    with utils.openw(args['rxl']) as out:
        n = len(fna.values()[0])
        out.write( str(len(fna))+" "+str(n)+"\n" )

        for k,v in fna.items():
            if len(k) > 14:
                k = k[:14]
            out.write( str(k)+" "*(15-len(str(k)[1:]))+str(v.seq) +"\n" )
