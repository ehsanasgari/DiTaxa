#!/usr/bin/env python

from __future__ import with_statement

import sys,argparse
import utils 
import collections

def read_params(args):
    parser = argparse.ArgumentParser(description='NCBI blast outfmt6 output processing')
    
    parser.add_argument('inp_f', metavar='INPUT_FILE', nargs='?', default=None, type=str, 
                        help="the input file [stdin if not present]")    
    parser.add_argument('out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str, 
                        help="the output file [stdout if not present]")    
    parser.add_argument('--pid', metavar='percentage identity', default=0, type=float,
                        help="minimum percentage identity (default 0)")
    parser.add_argument('--length', metavar='length', default=0, type=int,
                        help="minimum alignment length (default 0)")
    parser.add_argument('--evalue', metavar='evalue', default=float(sys.maxint), type=float,
                        help="maximum expect value (default INF)")
    parser.add_argument('--bitscore', metavar='bitscore', default=0.0, type=float,
                        help="minimum bit score (default 0)")
    parser.add_argument('--pid_col', metavar='pid_col', default=3, type=int,
                        help="column of the pid value (default 3)")
    parser.add_argument('--length_col', metavar='length_col', default=4, type=int,
                        help="column of the length value (default 4)")
    parser.add_argument('--evalue_col', metavar='evalue_col', default=11, type=int,
                        help="column of the expect value (default 11)")        
    parser.add_argument('--bitscore_col', metavar='bitscore_col', default=12, type=int,
                        help="column of the bit score (default 12)")
    parser.add_argument('-n', metavar='number of top hits', default=-1, type=int,
			help="number of results to show (default -1 means all)")
    parser.add_argument('-t', metavar='number of top hits per unique query', default=-1, type=int,
			help="number of results per unique qeury to show (default -1 means all)")
    parser.add_argument('-s', default=None, choices=["evalue","bitscore","length","pid"], type=str,
			help="the measure used to sort the output (default bitscore)")
    return vars(parser.parse_args()) 

def blast_ncbi_outfmt6_screen(par):
    finp,fout = bool(par['inp_f']), bool(par['out_f'])

    inp_mat = (l.rstrip('\n').split("\t") for l in (utils.openr(par['inp_f']) if finp else sys.stdin))

    out_mat =(l for l in inp_mat 
                    if float(l[par['pid_col']-1]) >= par['pid'] and
                       float(l[par['length_col']-1]) >= par['length'] and
                       float(l[par['evalue_col']-1]) <= par['evalue'] and
                       float(l[par['bitscore_col']-1]) >= par['bitscore']  )

    if 's' in par and par['s']:
        if par['s'] == 'pid':
            col = par['pid']-1
        elif par['s'] == 'evalue':
            col = par['evalue_col']
        elif par['s'] == 'length':
            col = par['length_col']-1
        elif par['s'] == 'bitscore':
            col = par['bitscore_col']-1

        out_mat = sorted( out_mat, 
                          key=lambda x: float(x[col-1]) )

        if 'n' in par and par['n'] > -1:
            out_mat = out_mat[:par['n']]
    
    unique_queries = collections.defaultdict( int ) 
    with utils.openw(par['out_f']) if fout else sys.stdout as out_file:
        if 't' in par and par['t'] > -1:
            for l in out_mat:
                unique_queries[l[0]] += 1
                if unique_queries[l[0]] > par['t']:
                    continue
                out_file.write("\t".join(l)+"\n")
        else:
            for l in out_mat:
                out_file.write("\t".join(l)+"\n")


if __name__ == '__main__':
    params = read_params(sys.argv)
    
    blast_ncbi_outfmt6_screen(params)
