#!/usr/bin/env python2
# -*- coding: utf8 -*-
## @file integration_plot_tutorial.py Advanced tutorial
# on how to examine the ensemble before and after
# the integration to plot meaningful results
# 
# refer to @ref TutorialPythonIntegrationPlot for the formatted version

# @page TutorialPythonIntegrationPlot Advanced Python Integration Tutorial
#
# In this tutorial, we explore the different ways to 
# examine an ensemble before and after integration.
# We also demonstrate how to convert the position and
# velocity of planets into orbital elements and finally,
# we plot the difference in semi-major axis of planets.
#
# Source code for this tutorial can be found at @ref py/integration_plot_tutorial.py
# \section setup Set-up
#
# The set-up is quite similar to approach taken in @ref TutorialPythonResume. 
#
# Import the packages and parse the command-line arguments

import swarmng
from os import path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config" , help="Config file", required = True)
args = parser.parse_args()

# Next, load the config file and the initial conditions
cfg = swarmng.Config.load(args.config)
 
fn = cfg["input"]
ext = path.splitext(fn)[1]
if ext == "txt" : ens = swarmng.DefaultEnsemble.load_from_text(fn)
else : ens = swarmng.DefaultEnsemble.load_from_bin(fn)

#
# Integration is the same as the basic tutorial
ref = ens.clone()
swarmng.init(cfg)
integ = swarmng.Integrator.create( cfg )
integ.ensemble = ens
integ.destination_time = float(cfg["destination_time"])
integ.integrate()

# \section ed Examining data
# 
# At this point, the integration is completed. For analysis
# we have two ensembles: `ref`, the reference ensemble
# that keeps a copy of initial conditions and `ens`, 
# the working copy that was simulated to the `destination_time`.
#
# In this example, we want to find the semi-major axis of planets and
# find the evolution of this variable over time.
# 
# The basic access methods only give us position and velocity of planets.
# To calculate the orbital elements, we use the function 
# @ref swarmng.keplerian_for_cartesian. 
#
# The function `extarct_semi_major` takes an ensemble
# and calculate the semi-major axis for all planets in
# all of systems and returns a two-dimentional array of
# the results. For easier plotting, the first dimension
# is the planet number and the second dimension is the 
# system number. 
#
# For large two-dimensional arrays, we use NumPy. It is
# much faster than regular Python arrays.
#
# The core function to calculate the semi-major axis is 
# @ref swarmng.keplerian_for_cartesian. It takes two 3-tuples
# for the information of the planet and the star and returns
# a tuple of all orbital elemets. For more details refer to 
# documentation for @ref swarmng.keplerian_for_cartesian.
import numpy

def extract_semi_major(en):
    shape = (ens.nbod-1,ens.nsys)
    value = numpy.empty(shape)
    for sysid, sys in enumerate(en):
        star = sys[0]
        center = (star.pos,star.vel,star.mass)

        for bodid in range(1,en.nbod):
            bod = sys[bodid]
            planet = (bod.pos,bod.vel,bod.mass)
            orbital_elements = swarmng.keplerian_for_cartesian(planet,center)
            value[bodid-1,sysid] = orbital_elements.a

    return value    

# We calculate semi-major axis values for both ensembles: initial and
# final conditions (final means after the simulation).
initial = extract_semi_major(ref)
final   = extract_semi_major(ens)
#
# The results are basically two-dimensional arrays of same shape.
# The next line normalizes the changes on the values. 
# Note
# that the operations are carried out element-wise, so that 
# the line of code here is equivalent to executing
# `change[i,j] = (final[i,j]-initial[i,j])/initial[i,j]` for 
# all valid `i`,`j`. 
change = (final-initial)/initial

#
# \section plt Plotting
#
# For every planet, we plot a scatter plot of change of semi-major axis
# where each data point represents the change in one of the systems.
# The plots from different planets are overlapped but distinguished 
# by different markers.
#
# X-axis represents the initial value of the semi-major axis calculated
# from the initial conditions.
#
# Y-axis represents the normalized change in the semi-major axis for 
# individual planets in different systems.
# 
markers = [ '.','x', '*', '^', 's', 'o' ]
import matplotlib.pyplot as P
for i in range(1,ens.nbod):
   P.scatter(initial[i-1],change[i-1],marker= markers[i-1], label="Planet {0}".format(i)) 

P.xlabel('Initial semi-major axis')
P.ylabel('Relative change of semi-major axis')
P.legend(title="Semi-major axis change")
P.savefig("semi_major_axis_comparison")
# 
# The generated plot is saved as `semi_major_axis_comparison` 
# in the currenty working directory.
#
# \section conc Conclusion
#
# In this, tutorial we explored how to extract information from an ensemble
# in a real application. You can change the tutorial to plot other
# orbital elements, or change the simulation part to more complicated one.
#
# We used NumPy and MatPlotLib packages, it is not possible to cover
# those packages in this tutorial. It is best to look for comprehensive
# tutorials of these packages on their respective websites.
# 
#
# **Where to go from here:**
# * @ref swarmng.keplerian_for_cartesian documentation for using other orbital element parameters.
# * @ref swarmng.Body documentation for other attributes of bodies that can be extracted from the ensemble
# * <a href="http://wiki.scipy.org/Tentative_NumPy_Tutorial">
#   Tentative NumPy Tutorial</a> for explanation of NumPy objects and
#   operators defined on them. There is also examples on element-wise 
#   operations on n-dimensional arrays.
# * <a href="http://matplotlib.org/users/pyplot_tutorial.html">
#   MatPlotLib tutorial</a> for customizing plots and creating other types of plots.
#  
