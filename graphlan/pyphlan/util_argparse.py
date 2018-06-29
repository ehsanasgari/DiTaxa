import argparse as ap
import sys


class util_argparse:

    def __init__( self, script_name = None ):
        self.script_name = script_name
        self.parser = ap.ArgumentParser()
        self.arg = self.parser.add_argument
        self.args = {}

    def inp_0( self ):
        self.arg( 'inp', metavar='INPUT_FILE', type=str, nargs='?', default=None )

    def out_0( self ):
        self.arg( 'out', metavar='OUTPUT_FILE', type=str, nargs='?', default=None )

    def version( self ):
        arg( '-v','--version', action='store_true', 
             help="Prints the current "+self.script_name+" version and exit\n" if self.script_name else "" ) 

    def get_inp_0( self ):
        if not 'inp' in self.args:
            return None
        if self.args['inp'] is None:
            return sys.stdin
        return self.args['inp']
    
    def get_out_0( self ):
        if not 'out' in self.args:
            return None
        if self.args['out'] is None:
            return sys.stdout
        return self.args['out']

    def parse( self ):
        self.args = vars(self.parser.parse_args())
    


