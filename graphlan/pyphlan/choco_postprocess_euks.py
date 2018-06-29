#!/usr/bin/env python

import sys
import collections
import utils
import pyphlan as ppa  

try:
    import argparse as ap
    import bz2 
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Profile ChocoPhlAn genes\n')

    p.add_argument( '--markers', required = True, default=None, type=str )
    p.add_argument( '--euk_taxonomy', required = True, default=None, type=str )
    #p.add_argument( '--mic_taxonomy', required = True, default=None, type=str )
    p.add_argument('out', nargs='?', default=None, type=str,
            help=   "the output summary file\n"
                    "[stdout if not present]")

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )

    cen2data = {}
    etree = ppa.PpaTree( args['euk_taxonomy'] )
    #vtree = ppa.PpaTree( args['vir_taxonomy'] )
    #mtree = ppa.PpaTree( args['mic_taxonomy'] )
    #vall = set([a.name for a in vtree.tree.get_terminals()])
    #mall = set([int(a.name[3:]) for a in mtree.tree.get_terminals()])
    eall = set([a.name for a in etree.tree.get_terminals()])

    lin = (l.split('\t') for l in open(args['markers']))
    for d1,cen,taxa,tax,d2,d3,coreness,d4,d5,d6,d7,tin,tout,d8,tsin,tsout in lin:
        tsin,tsout = set(["t__"+a for a in tsin.strip().split(":") if a]),set(["t__"+b for b in tsout.strip().split(":") if b])
        #tsin,tsout = set([int(a) for a in tsin.strip().split(":") if a]),set([int(b) for b in tsout.strip().split(":") if b])
        cen2data[int(cen)] = {'taxa':int(taxa),'tax':tax,'coreness':float(coreness),'tsin':tsin,'tsout':tsout}


    tax2cen = collections.defaultdict( set )
    for k,v in cen2data.items():
        tax2cen[v['tax']].add( k )
   
    for t,cs in tax2cen.items():
        q0 = sorted([c for c in cs if not cen2data[c]['tsout']],key=lambda x:-cen2data[x]['coreness'])

        """
        ok,err = [],[] 
        for c in cs:
            if c in q0:
                continue
            vext = vall & cen2data[c]['tsout']
            if vext: 
                err.append( c )
            else:
                ok.append( c )

        q10,q20 = {},{}
        for c in ok:
            q10[c] = mtree.lca( [str(tt) for tt in cen2data[c]['tsout'] if str(tt) in mall] ).full_name.split(".t__")[0]
        for c in err:
            q20[c] = mtree.lca( [str(tt) for tt in cen2data[c]['tsout'] if str(tt) in mall] ).full_name.split(".t__")[0]
        """

        def maxn(x,y): return x > 2000 

        n = 0
        for o in q0:
            sys.stdout.write( "\t".join([ str(o),str(t),""]) + "\n" )
            n += 1
            if n > 2000:
                break

        if maxn( n, len(cs) ):
            continue
        
        """
        added = set()
        for taxlev in "sgfoc":
            for c in ok:
                if True: # taxlev+"__" in q10[c]:
                    if c not in added:
                        sys.stdout.write( "\t".join([str(c),str(t),str(q10[c])]) + "\n"  )
                        added.add(c)
                        n += 1
                if maxn( n, len(cs) ):
                    break
            for c in err:
                if True: # taxlev+"__" in q20[c]:
                    if c not in added:
                        sys.stdout.write( "\t".join([str(c),str(t),str(q20[c])]) + "\n"  )
                        added.add(c)
                        n += 1
                if maxn( n, len(cs) ):
                    break
            if maxn( n, len(cs) ):
                break
        """
        ##print "+++++++++++++++++++++"

        ##print t,len(cs),len(q0),err
