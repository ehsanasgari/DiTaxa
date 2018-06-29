#!/usr/bin/env python

import util_argparse as ap
from tab import *

def read_params():
    p = ap.util_argparse() 
    p.inp_0()
    p.out_0()
    p.arg( '-r','--row',type=str, default = None)  
    p.arg( '-c','--col',type=str, default = None)  
    p.arg( '-v','--value',type=str, default = None)  
    p.arg( '-i','--inverse', action='store_true', default = False)  
    p.arg( '--regex',type=str, default = None)  
    p.arg( '--new_val',type=str, default = None)  
    p.arg( '-s','--sep',type=str, default = "\t")  
    p.parse()
    return p

if __name__ == '__main__':
    parser = read_params()
    args = parser.args
    
    tab = Tab( fin = parser.get_inp_0(), sep = args['sep'] ) 
    tab.sub( args['row'], args['col'], args['value'], args['regex'], args['new_val'], args['inverse'] )
    tab.save( parser.get_out_0() )

