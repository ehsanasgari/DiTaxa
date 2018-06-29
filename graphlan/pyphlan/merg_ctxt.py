#!/usr/bin/env python

import sys
import collections

try:
    import argparse as ap
    import bz2 
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description="Merge two core txt files"
            "where the second are subclusters of the first" )

    p.add_argument( 'ctxt', default=None, type=str,
            help=   "the main ctxt")
    p.add_argument('--out_ctxt', default=None, type=str,
            help=   "the output txt file compressef if given with bz2 extension")
    p.add_argument('sctxt', nargs='*', default=None, type=str,
            help=   "the list of subclustere ctxt")

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )

    openr = bz2.BZ2File if args['ctxt'].endswith(".bz2") else open
    mtxt = dict([([v[0],(v[1:] if len(v) > 1 else [])]) for v in (l.strip().split('\t') for l in openr(args['ctxt']))])
    
    for f in args['sctxt']:
        openr = bz2.BZ2File if f.endswith(".bz2") else open
        ctxt = dict([([v[0],(v[1:] if len(v) > 1 else [])]) for v in (l.strip().split('\t') for l in openr(f))])

        for k,v in mtxt.items():
            nv = v
            for kk in [k]+v:
                if kk in ctxt:
                    nv += ctxt[kk]
            mtxt[k] = nv

    openw = bz2.BZ2File if args['out_ctxt'].endswith(".bz2") else open
    with (openw(args['out_ctxt'],"w") if args['out_ctxt'] else sys.stdout) as out:
        for k,v in sorted(mtxt.items(),key=lambda x:-len(x[1])):
            out.write( "\t".join([k]+list(v)) +"\n" )
