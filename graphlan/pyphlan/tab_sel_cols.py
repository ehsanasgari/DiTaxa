#!/usr/bin/env python

import util_argparse as ap
from tab import *

def read_params():
    p = ap.util_argparse() 
    p.inp_0()
    p.out_0()
    p.arg( '-r','--row',type=str, default = None)  
    p.arg( '-v','--value',type=str, default = None)  
    p.arg( '--regex',type=str, default = None)  
    p.arg( '-s','--sep',type=str, default = "\t")  
    p.parse()
    return p

if __name__ == '__main__':
    parser = read_params()
    args = parser.args
    
    tab = Tab( fin = parser.get_inp_0(), sep = args['sep'] ) 
    tab.sel_columns( args['row'], args['value'], args['regex'] )
    tab.save( parser.get_out_0() )

