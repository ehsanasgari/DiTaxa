#!/usr/bin/env python

import sys
import collections
import utils
import pandas
import StringIO
import tempfile

try:
    import argparse as ap
    import bz2 
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Generate a taxonomic txt hierarchy'
            ' with genome assignments from IMG taxa table\n')

    
    p.add_argument( 'img', nargs='?', default=None, type=str,
            help=   "the input uc file [stdin if not present]")
    p.add_argument('txt', nargs='?', default=None, type=str,
            help=   "the output txt file compressed if fiven with bz2 extension\n"
                    "[stdout if not present]")
    p.add_argument('--NCBI_names', metavar="NCBI names.dmp", default = None )
    p.add_argument('--corrections', metavar="Correction file", default = None )
    p.add_argument('-d', metavar="Domain",
            default='Mic', choices=['Mic','Vir','Euk'] )

    return vars( p.parse_args() )

qm = "?"
def get_qm( s, ncbiid, t, tl, existing_species ):
    name = t['Genome Name / Sample Name'].replace("Candidatus ","").split(" ")[:2]
    genus, species = name if len(name) > 1 else ("","")
    if s is None or s == "-1" or not s or type(s) != str or  s in ['Unclassified','unclassified'] or s.split("_")[0] in ['bacterium','sp','sp.','Sp','Sp.','spp','spp.','Spp','Spp.']:
        if tl == 's' and species in existing_species[t['Genus']] and genus == t['Genus']:
            s = species
            #print str(t.name)+"\t+\t"+str(genus)+" "+str(species)
        else: # species not in ['sp','sp.','Sp','Sp.','spp','spp.','Spp','Spp.']:
            #print "--",species
            if tl == 's' and genus == t['Genus']  and genus + " " + species == str(ncbiid[t['NCBI Taxon ID']]) and species not in ['bacterium','sp','sp.','Sp','Sp.','spp','spp.','Spp','Spp.']:
                s = species
                #print str(t.name)+"\t**\t"+str(genus)+" "+str(species)
            else:
                #print ">", genus, t['Genus'],  t['NCBI Taxon ID'], genus + " " + species, str(ncbiid[t['NCBI Taxon ID']])
                return qm
    return s.replace("Candidatus ","").replace(" ","_").replace(".","_").replace(",","_")

def add_unnamed( arr ):
    gunamed = "?" in arr[-2]
    if "?" in arr[0]:
        return a
    lt,last = 'd',arr[0]
    v = [arr[0]]
    for a in arr[1:-2]:
        if "?" in a and not gunamed:
            v.append( a[:3] + last+"_noname" )
        else:
            v.append( a )
        lt,last = a[0],(a[3:] if "?" not in a else last)
    v += arr[-2:] 
        
    return v 

if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )

    tax_lev = "dpcofgs"
    tax_lev_exp = ['Domain','Phylum','Class','Order','Family','Genus','Species']

    fp = tempfile.TemporaryFile()

    if args['corrections']:
        with utils.openr(args['corrections']) as inp:
            frto = {}
            frtoid = {}
            for pat in (l.split('\t') for l in inp):
                if len(pat) == 2:
                    frto[pat[0].strip()] = pat[1].strip()
                else:
                    frtoid[pat[0].strip()] = (pat[1].strip(),pat[2].strip())
            with utils.openr(args['img'],"rU") as inpf:
                nfa = []
                for l in inpf:
                    nf = l
                    for f,t in frto.items():
                        nf = nf.replace(f,t)
                    for i,(f,t) in frtoid.items():
                        if l.startswith(i+"\t"):
                            nf = nf.replace(f,t)
                    nfa.append(nf) 
                fp.write( "".join( nfa ) )
                fp.seek(0)

    with (utils.openr(args['img'],"rU") if not args['corrections'] else fp) as inp:
        table = pandas.read_table( inp, sep='\t', index_col=0)
 
    if args['d'] == 'Mic': 
        table = table[table['Domain'].isin(['Bacteria','Archaea'])]
        table = table[table['Gene Count'] > 250.0]
        table = table[table['Genome Size'] > 50000.0]
    elif args['d'] == 'Vir':
        table = table[table['Domain'].isin(['Viruses'])]
        table = table[table['Gene Count'] > 0.0]
    elif args['d'] == 'Euk':
        table = table[table['Domain'].isin(['Eukaryota'])]
        able = table[table['Gene Count'] > 1000.0]
        table = table[table['Genome Size'] > 500000.0]

    toexcl = ['Candidatus','sp','sp.','Sp','Sp.','spp','spp.','Spp','Spp.','Unclassified','unclassified',"-1"]
    #table['Species'] = [("-1" if t in toexcl else t) for t in table['Species']]

    infos = ['Genome Name / Sample Name','NCBI Taxon ID']
    table = table.reindex(columns=infos+tax_lev_exp+['Genome Name'])
    table['Genus'] = [t.replace("Candidatus ","") for t in table['Genus']]

    existing_species_l = [(b,a) for a,b in list(set(zip(list(table['Species']),list(table['Genus'])))) if a not in toexcl]
    existing_species = collections.defaultdict( set )
    for a,b in existing_species_l:
        existing_species[a].add(" ".join(b.replace("Candidatus ","").split(" ")[:2]))

    ncbi_species = set()
    ncbiid = dict([(int(a),None) for a in table['NCBI Taxon ID']])

    if args['NCBI_names']:
        for line in (l.strip().split("|") for l in open(args['NCBI_names'])):
            if int(line[0]) in ncbiid and "scientific name" in line[3]:
                ncbiid[int(line[0])] = " ".join(line[1].replace("Candidatus ","").strip().split(" ")[:2])

    with utils.openw(args['txt']) as out:
        for i,t in table.iterrows():
            out.write( "\t".join( [ str(-int(i)),
                                  ".".join( add_unnamed(["__".join([taxl, get_qm(t[taxle],ncbiid,t,taxl,existing_species)]) 
                                      for taxl,taxle in  zip(tax_lev,tax_lev_exp)]) + ["t__"+str(-int(i))]
                                        ),
                                  #str(t['Genome Name'])
                                  ] ) + "\n")
    
    
