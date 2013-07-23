# @file stats_tutorial.py Tutorial on extracting statistical information
# from a data file

# @page TutorialPythonStats Extracting statistical information from data files
#
# First try to make sure that correct command-line arguments
# are provided.
import swarmng
from sys import argv,exit
from os import path
if len(argv) != 2 or not path.exists(argv[1]):
    print("Usage: {0} <path_to_data_file>".format(argv[0]))
    exit(1)

# Next load the data file based on the extension
fn = argv[1]
ext = path.splitext(fn)[1]
if ext == ".txt" : ens = swarmng.DefaultEnsemble.load_from_text(fn)
else : ens = swarmng.DefaultEnsemble.load_from_bin(fn)
#
#
# Now gather all data in a 3-dimensional NumPy array.
# NumPy has internal statstical information analysis
# functions that makes it easier to find different
# statistical information.
#
# Dimensions are:
# * planet number 0 to nbod-1
# * orbital element number : 0 to 5
# * system id
#
import numpy

shape = (ens.nbod-1, 6, ens.nsys)
value = numpy.empty(shape)
for sysid, sys in enumerate(ens):
    star = sys[0]
    center = (star.pos,star.vel,star.mass)

    for bodid in range(1,ens.nbod):
        bod = sys[bodid]
        planet = (bod.pos,bod.vel,bod.mass)
        orbital_elements = swarmng.keplerian_for_cartesian(planet,center)
        value[bodid-1, :,sysid] = orbital_elements

# Now iterate through the 3-dimensional array
# and calculate the median of orbital element
# values for each body.
#
for bodid in range(1,ens.nbod):
    orbital_elements = [0,0,0,0,0]
    for i in range(0,5):
        orbital_elements[i] = numpy.median(value[bodid-1,i,:])
    print("Body {0}: {1}, {2}, {3}, {4}, {5}".format(bodid,*orbital_elements))


