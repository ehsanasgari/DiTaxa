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
    p = ap.ArgumentParser(description='Convert core gene files to core gene summaries\n')

    p.add_argument( 'cg', nargs='?', default=None, type=str,
            help=   "the input cg file [stdin if not present]")
    p.add_argument('cgs', nargs='?', default=None, type=str,
            help=   "the output summary file\n"
                    "[stdout if not present]")

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )

    gid2cores = collections.defaultdict( set )
    #with (open(args['uc']) if args['uc'] else sys.stdin) as inp:
    with utils.openr(args['cg']) as inp:
        for line in (l.split('\t') for l in inp):
            if int(line[0]) > 0:
                gid,clade,ncore,ngenomes,pv =  line[:5]
            else:
                gid,clade,ncore,ngenomes,pv =  line[1:6]
            gid2cores[gid].add( (clade,ncore,ngenomes,pv) )

    clades2cores = collections.defaultdict( set )
    for k,v in gid2cores.items():
        if len(v) > 1:
            continue
        clades2cores[list(v)[0][0]].add( k )

    #openw = bz2.BZ2File if args['txt'].endswith(".bz2") else open
    with utils.openw(args['cgs']) as out:
        for k,v in sorted(clades2cores.items(),key=lambda x:-len(x[1])):
            out.write( "\t".join([k,str(len(v))]) +"\n" )

