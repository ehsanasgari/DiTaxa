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
    p.add_argument('-n', metavar="Minimum number of matching taxa",
            default=0, type=int )
    p.add_argument('-m', metavar="Minimum number of matching taxa",
            default=0, type=int )
    p.add_argument('-p', metavar="Prefix for taxon names",
            default="", type=str )
    p.add_argument('--sk', action='store_true' )
    p.add_argument('-s', action='store_true' )
    p.add_argument('--original_gene_names', action='store_true' )

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )
    
    gint = str if args['sk'] else int

    valin, kin = [], []
    with utils.openr( args['ctxt'] ) as inp:
        for l in inp:
            tset = set([gint(a) for a in l.strip().split('\t')][1:])
            if len(tset) < args['n']:
                continue
            valin.append(tset)
            kin.append( l.strip().split('\t')[0] )
    all_t = set()
    for v in valin:
        all_t |= v

    res = {}
    for t in all_t:
        #if len(t) < args['n']:
        #    continue
        res[t] = [int(t in v) for v in valin]

    if args['s']:
        keys = sorted( res.keys() )
        v = ["".join([str(vvv) for vvv in vv]) for vv in zip(*[res[k] for k in keys])]
        resk = collections.defaultdict( int )
        for vv in v:
            resk[ vv ] += 1
        res = [r for r in [ [int(kk)*v for kk in k] for k,v in resk.items()] if max(r) > args['m']]
        res = dict([(keys[i],v) for i,v in enumerate(zip(*res))])

    with utils.openw(args['txt']) as out:
        n = len(res.values()[0])
        n_s = int(float(n)*args['subsample'])
        indok = set(random.sample( list(range(n)), n_s))

        torem = set()
        if not args['s'] and args['m']:
            torem = set([i for i,l in enumerate(zip(*res.values())) if sum(l) < args['m']])

        if args['original_gene_names']:
            out.write( "\t".join(['genes']+[kin[v] for v in range(n)]) + "\n" )
        else:
            out.write( "\t".join(['genes']+["g"+str(v) for v in range(n)]) + "\n" )

        for k,v in res.items():
            out.write( args['p'] + str(k)+"\t"+"\t".join([str(s) for i,s in enumerate(v) if i in indok and i not in torem]) +"\n" )









