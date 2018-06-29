from pandas import *
import sys
import re
import msg

class Tab:

    def __init__( self, fin = sys.stdin, sep = "\t" ):
        self.data = read_csv( fin, sep=sep, index_col=0 ) 

    def sel_rows( self, col = None, val = None ):
        #self.data = self.data[self.data[col] == val]
        pass     

    def eq( self, x, y ):
        t = type(y)
        if x == y:
            return True
        try:
            return t(x) == y
        except Exception:
            return False

    def get_cols( self, row = None, val = None, regex = None ):
        cols = []
        if row is None:
            if regex:
                for c,cs in self.data.iteritems():
                    if any((re.search(regex,str(v)) for v in cs)):
                        #del self.data[c]
                        cols.append(c)
            elif val:
                for c,cs in self.data.iteritems():
                    if any((val == str(v) for v in cs)):
                    #if val not in cs:
                        #del self.data[c]
                        cols.append(c)
            else:
                msg.exit("Error") 
        else:
            if regex:
                for c,cs in self.data.iteritems():
                    if re.search(regex,str(cs[row])):
                        #del self.data[c]
                        cols.append(c)
            elif val:
                for c,cs in self.data.iteritems():
                    if self.eq(val, cs[row]):
                        #del self.data[c]
                        cols.append(c)
            else:
                msg.exit("Error") 
        return cols

    def sel_columns( self, row = None, val = None, regex = None ):
        cols = set(self.get_cols( row, val, regex ))
        for c in set(self.data.columns) - cols:
            del self.data[c]

    def sub( self, row = None, col = None, val = None, regex = None, new_val = None, inverse = False ):
        if row and col:
            if regex:
                if inverse:
                    if not re.match( regex, str(self.data[col][row]) ):
                        self.data[col][row] = new_val
                else:
                    self.data[col][row] = re.sub( regex, new_val, str(self.data[col][row]) )
            else:
                self.data[col][row] = new_val
        elif col:
            if regex:
                if inverse:
                    for i,r in enumerate(self.data[col]):
                        if not re.match( regex, str(r) ):
                            self.data[col][i] = new_val
                else:
                    for i,r in enumerate(self.data[col]):
                        self.data[col][i] = re.sub( regex, new_val, str(r) )
            elif val:
                for i,r in enumerate(self.data[col]):
                    if str(r) == val:
                        self.data[col][i] = new_val
            else:
                msg.exit("Error")
        elif row:
            if regex:
                if inverse:
                    for i,c in enumerate(self.data.ix[row]):
                        if not re.match( regex, str(c) ):
                            self.data.ix[row][i] = new_val
                else:
                    for i,c in enumerate(self.data.ix[row]):
                        self.data.ix[row][i] = re.sub( regex, new_val, str(c) )
            elif val:
                for i,c in enumerate(self.data.ix[row]):
                    if str(c) == val:
                        self.data.ix[row][i] = new_val
            else:
                msg.exit("Error")
        
        #cols = self.get_cols( row, val, regex )
         

    def save( self, outf = sys.stdout, sep = "\t" ):
        self.data.to_csv(outf,sep="\t")

