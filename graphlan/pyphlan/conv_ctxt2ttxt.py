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

    p.add_argument( 'ctxt', nargs='?', default=None, type=str,
            help=   "the input ctxt file [stdin if not present]")
    p.add_argument('--g2t', metavar="Mapping file from genes to taxa",
            default=None, type=str )
    p.add_argument('-n', metavar="Minimum size of cluster to report (def 1, means all)",
            default=1, type=int )
    p.add_argument('--t2g', metavar="Mapping file from taxa to genes",
            default=None, type=str )
    p.add_argument('--sk', action='store_true' )
    p.add_argument('txt', nargs='?', default=None, type=str,
            help=   "the output gtxt file, compressed if fiven with bz2 extension\n"
                    "[stdout if not present]")

    return vars( p.parse_args() )


if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )

    gint = str if args['sk'] else int

    """
    with utils.openr( args['ctxt'] ) as inp:
        valin = [[int(a) for a in l.strip().split('\t')] for l in inp]
    if not args['g2t'] and not args['t2g']:
        sys.stdout.write("Error one of --t2g and --g2t must be provided\n")
        sys.exit(0)
    g2t = {}


    if args['g2t']:
        with utils.openr( args['g2t'] ) as inp:
            #g2t = dict(([int(a) for a in l.strip().split('\t')] for l in inp))
            for l in inp:
                f,t = l.strip().split('\t')
                g2t[int(f)] = int(t)
    elif args['t2g']:
        with utils.openr( args['t2g'] ) as inp:
            for ll in (l.strip().split('\t') for l in inp):
                for g in ll[1:]:
                    g2t[int(g)] = int(ll[0])
    
    for v in valin:
        valin[j] = [v[0]]+list(set(g2t[vv] for vv in v))

    with utils.openw(args['txt']) as out:
        for v in sorted(valin,key=lambda x:-len(x)):
            out.write( "\t".join([str(s) for s in v]) +"\n" )

    """

    if not args['g2t'] and not args['t2g']:
        sys.stdout.write("Error one of --t2g and --g2t must be provided\n")
        sys.exit(0)
    g2t = {}


    if args['g2t']:
        with utils.openr( args['g2t'] ) as inp:
            #g2t = dict(([int(a) for a in l.strip().split('\t')] for l in inp))
            for l in inp:
                f,t = l.strip().split('\t')
                g2t[gint(f)] = gint(t)
    elif args['t2g']:
        with utils.openr( args['t2g'] ) as inp:
            for ll in (l.strip().split('\t') for l in inp):
                for g in ll[1:]:
                    g2t[gint(g)] = gint(ll[0])
    
    #for v in valin:
    #    valin[j] = [v[0]]+list(set(g2t[vv] for vv in v))

    with utils.openw(args['txt']) as out:
        with utils.openr( args['ctxt'] ) as inp:
            for l in inp:
                valin = [gint(a) for a in l.strip().split('\t')]
                valin = [valin[0]]+list(set(g2t[vv] for vv in valin)) 
    
                #for v in sorted(valin,key=lambda x:-len(x)):
                if len(valin) > args['n']:
                    out.write( "\t".join([str(s) for s in valin]) +"\n" )



