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
    p.add_argument( '--vir_taxonomy', required = True, default=None, type=str )
    p.add_argument( '--mic_taxonomy', required = True, default=None, type=str )
    p.add_argument('out', nargs='?', default=None, type=str,
            help=   "the output summary file\n"
                    "[stdout if not present]")

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )

    cen2data = {}
    vtree = ppa.PpaTree( args['vir_taxonomy'] )
    mtree = ppa.PpaTree( args['mic_taxonomy'] )
    vall = set([a.name for a in vtree.tree.get_terminals()])
    #mall = set([int(a.name[3:]) for a in mtree.tree.get_terminals()])
    mall = set([a.name for a in mtree.tree.get_terminals()])

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

        ok,err = [],[] 
        for c in cs:
            if c in q0:
                continue
            #print "======================",list(vall)[:3], list(cen2data[c]['tsout'])[:3]
            vext = vall & cen2data[c]['tsout']
            if vext: 
                #print "ERRRRRRRRRRRRRRRRRRRR"
                err.append( c )
            else:
                #print "OOOOOOOOOOOOOOOOOR"
                ok.append( c )

        q10,q20 = {},{}
        for c in ok:

            q10[c] = mtree.lca( [str(tt) for tt in cen2data[c]['tsout'] if str(tt) in mall] ).full_name.split(".t__")[0]
            #print ["t__"+str(tt) for tt in cen2data[c]['tsout']]
            #print list(vall)[:3]
            #print ["t__"+str(tt) for tt in cen2data[c]['tsout'] if "t__"+str(tt) in mall]
            #print "q1000000000000",q10[c], [str(tt) for tt in cen2data[c]['tsout']]
        for c in err:
            q20[c] = mtree.lca( [str(tt) for tt in cen2data[c]['tsout'] if str(tt) in mall] ).full_name.split(".t__")[0]
            #print ["t__"+str(tt) for tt in cen2data[c]['tsout']]
            #print "q2000000000000", q20[c], ["t__"+str(tt) for tt in cen2data[c]['tsout']]

        def maxn(x,y): return x > 50 or ( x > 10 and y > 15 )

        #print q0,q10,q20

        n = 0
        for o in q0:
            sys.stdout.write( "\t".join([ str(o),str(t),""]) + "\n" )
            n += 1
            if n > 50:
                break

        if maxn( n, len(cs) ):
            continue
        ##print "============" 
        added = set()
        for taxlev in "sgfoc":
            ##print "a"
            for c in ok:
                ##print "q",taxlev+"__", q10[c] 
                if True: # taxlev+"__" in q10[c]:
                    ##print taxlev+"__" in q10[c]
                    if c not in added:
                        sys.stdout.write( "\t".join([str(c),str(t),str(q10[c])]) + "\n"  )
                        added.add(c)
                        n += 1
                if maxn( n, len(cs) ):
                    break
            ##print "b"
            for c in err:
                #print c,q20[c]
                if True: # taxlev+"__" in q20[c]:
                    if c not in added:
                        sys.stdout.write( "\t".join([str(c),str(t),str(q20[c])]) + "\n"  )
                        added.add(c)
                        n += 1
                if maxn( n, len(cs) ):
                    break
            if maxn( n, len(cs) ):
                break

        ##print "+++++++++++++++++++++"

        ##print t,len(cs),len(q0),err


    """
    cscores = collections.defaultdict( set )
    fwmarkers = {} 
    maps = collections.defaultdict( set )
    tree = ppa.PpaTree( args['taxonomy'] )
    clades2terms = ppa.clades2terms( tree.tree ) 

    clades2taxa = dict([(clade.full_name,set([-int(taxon.name[4:]) for taxon in taxa])) for clade,taxa in clades2terms.items()])
    ntaxa = 0
    for v in clades2terms.values():
        ntaxa += len(v)

    for l in open( args['cscores'] ):
        gene_seed, clade, n, n_tot, coreness = l.strip().split('\t')
        gene_seed, clade, n, n_tot, coreness = int(gene_seed), clade, int(n), int(n_tot), float(coreness)
        cscores[gene_seed].add( (clade, n, n_tot, coreness) )
    
    for i,l in enumerate(open( args['fwmarkers'] )):
        taxa_id, gene_seed, clade, n, n_tot, coreness, n_ext_seeds, n_ext_taxa, uniqueness, ext_taxa = l.strip().split('\t')
        taxa_id, gene_seed, n, n_tot, coreness, n_ext_seeds, n_ext_taxa, uniqueness = int(taxa_id), int(gene_seed), int(n), int(n_tot), float(coreness), int(n_ext_seeds), int(n_ext_taxa), float(uniqueness)
        fwmarkers[gene_seed] = (taxa_id, clade, n, n_tot, coreness, n_ext_seeds, n_ext_taxa, uniqueness)

    for l in open( args['maps'] ):
        line = l.strip().split('\t')
        if len(line) < 2:
            continue
        fr,to = line
        fr,to = int(fr.split("_")[0]), to # can be improved!!!!
        maps[fr].add( to )

    g2t,g2c,c2g = {},{},{}
    if args['g2t']:
        with open( args['g2t'] ) as inp:
            g2t = dict(([int(a) for a in l.strip().split('\t')] for l in inp))
    elif args['t2g']:
        with open( args['t2g'] ) as inp:
            for ll in (l.strip().split('\t') for l in inp):
                for g in ll[1:]:
                    g2t[int(g)] = int(ll[0])
        
    genomes = set([g2t[g] for g in cscores]) 

    with open( args['g2c'] ) as inp:
        for l in inp:
            line = list(l.strip().split('\t'))
            #if int(line[0]) not in genomes:
            #    continue
            #vals = [int(a) for a in line if utils.is_number(a)]
            vals = [a for a in line]
            if len(vals) > 1:
                g2c[int(vals[0])] = vals[1:]
    for g,c in g2c.items():
        for cc in c:
            c2g[cc] = g
    with utils.openw( args['out'] ) as out:
        for gene_seed,cscores_t in cscores.items():
            taxa = g2t[gene_seed]
            for clade, n, n_tot, coreness in cscores_t:
                out.write( "\t".join(["CSCORE",str(gene_seed),str(taxa),clade,str(n), str(n_tot), str(coreness)]) +"\n" )
            
            # anche sotto ???

            if gene_seed in fwmarkers:
                taxa_id, clade, n, n_tot, coreness, n_ext_seeds, n_ext_taxa, uniqueness = fwmarkers[gene_seed]
                if uniqueness < 0.01:
                     out.write( "\t".join(["FWMARKER",str(gene_seed),str(taxa),clade,str(n), str(n_tot), str(coreness),
                                                      str(n_ext_seeds), str(n_ext_taxa), str(1.0-uniqueness)]) +"\n" ) 

                if gene_seed in maps:
                    ext_tax = set([(c2g[s] if s in c2g else 0) for s in maps[gene_seed]])
                    ext_tax_ok = ext_tax & clades2taxa[clade]
                    ext_tax_ko = ext_tax - clades2taxa[clade]
                    ext_tax_okl = ":".join([str(s) for s in ext_tax_ok])
                    ext_tax_kol = ":".join([str(s) for s in ext_tax_ko])
                    muniq = 1.0 - float(len(ext_tax_kol)) / float(ntaxa)

                    out.write( "\t".join(["MARKER",str(gene_seed),str(taxa),clade,str(n), str(n_tot), str(coreness),
                                                   str(n_ext_seeds), str(n_ext_taxa), str(1.0-uniqueness), 
                                                   str(len(ext_tax)), str(len(ext_tax_ok)), str(len(ext_tax_ko)), 
                                                   str(muniq), ext_tax_okl, ext_tax_kol]) +"\n" ) 

    """
