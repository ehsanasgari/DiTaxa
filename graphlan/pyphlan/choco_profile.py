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

    p.add_argument( '--cscores', required = True, default=None, type=str )
    p.add_argument( '--fwmarkers', required = True, default=None, type=str )
    p.add_argument( '--maps', required = True, default=None, type=str )
    p.add_argument( '--taxonomy', required = True, default=None, type=str )
    p.add_argument('--g2t', metavar="Mapping file from genes to taxa", default=None, type=str )
    p.add_argument('--t2g', metavar="Mapping file from taxa to genes", default=None, type=str )
    p.add_argument('--g2c', metavar="Mapping file from genomes to contigs", required = True, default=None, type=str )
    p.add_argument('out', nargs='?', default=None, type=str,
            help=   "the output summary file\n"
                    "[stdout if not present]")

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )

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
    gid2cores = collections.defaultdict( set )
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

    with utils.openw(args['cgs']) as out:
        for k,v in sorted(clades2cores.items(),key=lambda x:-len(x[1])):
            out.write( "\t".join([k,str(len(v))]) +"\n" )
    """
