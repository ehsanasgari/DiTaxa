#!/usr/bin/env python

import sys
import collections
import utils
try:
    import argparse as ap
    import bz2 
    import random
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Convert Usearch ".uc" files in tab-delimited'
            ' files with the seed as first field followed by the other IDs\n')

    p.add_argument( 'ctxt', nargs='?', default=None, type=str,
            help=   "the input uc file [stdin if not present]")
    p.add_argument('txt', nargs='?', default=None, type=str,
            help=   "the output txt file compresse if fiven with bz2 extension\n"
                    "[stdout if not present]")
    p.add_argument('--subsample', metavar="Subsampling rate",
            default=1.0, type=float )
    p.add_argument('-n', action='store_true' )
    p.add_argument('-p', metavar="Prefix for taxon names",
            default="", type=str )
    p.add_argument('--sk', action='store_true' )

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )
    
    gint = str if args['sk'] else int

    valin = []
    with utils.openr( args['ctxt'] ) as inp:
        for l in inp:
            tset = set([gint(a) for a in l.strip().split('\t')][1:])
            valin.append(tset)
    all_t = set()
    for v in valin:
        all_t |= v

    res = {}
    for t in all_t:
        #if len(t) < args['n']:
        #    continue
        res[t] = [int(t in v) for v in valin]

    with utils.openw(args['txt']) as out:
        mat_res = {}
        for k1 in res:
            mat_res[k1] = {}
            for k2 in res:
                mat_res[k1][k2] = float(len([1 for x,y in zip(res[k1],res[k2]) if x == 1 and y == 1]))/max(sum(res[k1]),sum(res[k2]))

        keys = sorted(res)
        out.write( "\t".join( ["x"]+keys ) +"\n" )
        for k in keys:
            out.write( "\t".join( [k] + [str(mat_res[k][k2]) for k2 in keys] ) +"\n" )


        """
        n = len(res.values()[0])
        n_s = int(float(n)*args['subsample'])
        indok = set(random.sample( list(range(n)), n_s))
        out.write( "\t".join(['genes']+["g"+str(v) for v in range(n)]) + "\n" )

        for k,v in res.items():
            out.write( args['p'] + str(k)+"\t"+"\t".join([str(s) for i,s in enumerate(v) if i in indok]) +"\n" )
        """
