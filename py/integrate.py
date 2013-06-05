#!/usr/bin/env python
import sys, getopt
from sys import path

path.append('../lib')
from libswarmng_ext import *

# Integrate ensembles based on the config file and input ensemble
 
def main(argv):
    inputfile = ''
    outputfile = ''
    cfgfile = ''

    try:
        opts, args = getopt.getopt(argv,"hc:i:o:",["cfgfile","ifile=","ofile="])
    except getopt.GetoptError:
        print 'integrate.py -c <cfgfile> -i <inputfile> -o <outputfile>'
        sys.exit(2)
        
    for opt, arg in opts:
        if opt == '-h':
            print 'integrate.py -c <cfgfile> -i <inputfile> -o <outputfile>'
            sys.exit()
	elif opt in ("-c", "--cfgfile"):
	    cfgfile = arg
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    cfg = load_config( cfgfile )

    ref = generate_ensemble( cfg )
    ens = ref.clone()
    integ = Integrator.create( cfg )
    integ.ensemble = ens
    integ.destination_time = 1.0

    sync

    integ.integrate()
    ref.save_to_text( outputfile )

    sync    
    

if __name__ == "__main__":
   main(sys.argv[1:])
