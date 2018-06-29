import sys

def info( s ):
    sys.stdout.write( s )
    sys.stdout.flush()

def exit( s ):
    sys.stderr.write(s+"\n")
    sys.stderr.write("Exiting ... \n")
    sys.exit(1)


