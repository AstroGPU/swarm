# -*- coding: utf-8 *-*
## @file stats_tutorial.py Tutorial on extracting statistical information
# from a data file

# @page TutorialPythonStats Extracting statistical information from data files
#
# In this tutorial, we revisit loading ensembles and extracting
# information from an ensemble; however, this time we are
# building a small utility to view statistical information about
# the contents of the ensemble.
#
# Source code for this tutorial can be found at @ref py/stats_tutorial.py
#
# @section ai Arguments
# First try to make sure that correct command-line arguments
# are provided.
import swarmng
from sys import argv,exit
from os import path
if len(argv) != 2 or not path.exists(argv[1]):
    print("Usage: {0} <path_to_data_file>".format(argv[0]))
    exit(1)
    
# @section load Loading the ensemble
# Next load the data file based on the extension
fn = argv[1]
ext = path.splitext(fn)[1]
if ext == ".txt" : ens = swarmng.DefaultEnsemble.load_from_text(fn)
else : ens = swarmng.DefaultEnsemble.load_from_bin(fn)
#
# @section cc Compute the orbital elements
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

# @section pr Write the information to standard output
# Now iterate through the 3-dimensional array
# and calculate the average and standard deviation
# of each orbital element value for each body, then
# print it to the output.
#
print("Body number, Semi-major axis, Eccentricity, Inclination, Longtitude of the ascending node, Argument of periapsis, Mean anomaly")
for bodid in range(1,ens.nbod):
    orbital_elements = [0,0,0,0,0]
    o = "{0}".format(bodid)
    for i in range(0,6):
      o += ", {0}±{1}".format(numpy.average(value[bodid-1,i,:]),numpy.std(value[bodid-1,i,:]))
    print(o)
#
# @section conc Conclusion
# This concludes this tutorial, to learn more about the
# orbital elements and the calculation method 
# look at @ref swarmng.keplerian_for_cartesian

