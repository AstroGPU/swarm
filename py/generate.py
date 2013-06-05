#!/usr/bin/env python

import sys, getopt
from sys import path

path.append('../lib')
from libswarmng_ext import *

# A example for generating a ensemble to a text output
# from a config file
# Usage: generate.py -c <configfile> -o <outputfile>

def main(argv):
    configfile = ''
    outputfile = ''
    
    try:
        opts, args = getopt.getopt(argv,"hc:o:",["cfile=","ofile="])
    except getopt.GetoptError:
        print 'generate.py -c <configfile> -o <outputfile>'
        sys.exit(2)
        
    for opt, arg in opts:
        if opt == '-h':
            print 'generate.py -c <configfile> -o <outputfile>'
            sys.exit()
        elif opt in ("-c", "--cfile"):
            configfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    cfg = load_config( configfile )

    ref = generate_ensemble( cfg )
    ref.save_to_text( outputfile )

if __name__ == "__main__":
   main(sys.argv[1:])

