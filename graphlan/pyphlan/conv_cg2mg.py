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
    p = ap.ArgumentParser(description='Add markerness'
            ' to core gene file\n')

    p.add_argument( 'ctxt', nargs='?', default=None, type=str,
            help=   "the input ctxt file [stdin if not present]")
    p.add_argument('--b6o', metavar="The outfmt6 file for the cores",
            default=None, type=str )
    p.add_argument('-n', metavar="Total number of target sets (total targets in the b6o file if unspecified)",
            default=None, required = True, type=int )
    p.add_argument('--g2t', metavar="Mapping file from genes to taxa",
            default=None, type=str )
    p.add_argument('--t2g', metavar="Mapping file from taxa to genes",
            default=None, type=str )
    p.add_argument('mtxt', nargs='?', default=None, type=str,
            help=   "the output mtxt file, compressed if fiven with bz2 extension\n"
                    "[stdout if not present]")

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )
   
    if not args['g2t'] and not args['t2g']:
        sys.stdout.write("Error one of --t2g and --g2t must be provided\n")
        sys.exit(0)
    g2t = {}
    if args['g2t']:
        with open( args['g2t'] ) as inp:
            g2t = dict(([int(a) for a in l.strip().split('\t')] for l in inp))
    elif args['t2g']:
        with open( args['t2g'] ) as inp:
            for ll in (l.strip().split('\t') for l in inp):
                for g in ll[1:]:
                    g2t[int(g)] = int(ll[0])
    
    with utils.openr( args['ctxt'] ) as inp:
        valin = (l.strip().split('\t') for l in inp)

        g2c = collections.defaultdict( set )
        
        if args['b6o']:
            inp_mat = ((int(a),int(b)) for a,b in (l.rstrip('\n').split("\t")[:2] for l in utils.openr(args['b6o'])))
    
            #all_targets = set()
            for fr,to in inp_mat:
                #all_targets.add( to )
                if fr != to:
                    g2c[fr].add( to )

        n = args['n'] # if args['n'] else len(all_targets)
        n = float(n)
    
        with utils.openw(args['mtxt']) as out:
            last,lastv = "",[]
            outbuf = []
            gt = None
            for v in valin:
                gt = int(v[0])
                if last == gt:
                    lastv = ""
                    continue
                if lastv:
                    outbuf.append( lastv )
                last = gt
                lastv = v
            if last and last != gt:
                outbuf.append( lastv )
            for v in outbuf:
                fr = int(v[0])
                frt = g2t[fr]
                nu = len(g2c[fr]) if args['b6o'] else 0
                targett = list(set([g2t[s] for s in g2c[fr]])) if args['b6o'] else []
                targets = ":".join([str(s) for s in targett]) if targett else "-"
                uniqueness = round(float(len(targett)) / n,3)
                out.write( "\t".join([str(g2t[fr])]+v+[str(nu),str(len(targett)),str(uniqueness),targets]) +"\n" )



