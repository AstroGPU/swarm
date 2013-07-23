#!/usr/bin/env python2
## @file resume_tutorial.py Tutorial on resuming integration from saved snapshots
# 
# refer to @ref TutorialPythonResume for formatted version.


# @page TutorialPythonResume Resume an integration
# In this tutorial, we explore more details about data files and integration.
# 
# The goals of this utility is to:
# * Run a basic integration from a configuration file
# * Automatically save the output 
# * Automatically resume from the last saved output
# * Integrate for a specific period of time
#
# Detail:
# * Do the integration based on a cfg file, input file is required
# * Save the state of data file in a snapshot file
# * Find the median of integrated time from the ensemble and
# add the predefined amount to it
#
# Source code for this tutorial can be found at @ref py/resume_tutorial.py
# \section ai Arguments
#
# First, we parse the comand line arguments using argparse
#
# 
import swarmng, argparse
from os import path

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config" , help="Config file", required = True)
parser.add_argument("-d", "--duration"   , help ="Duration of integration in AU", default=10.0, type=float )
args = parser.parse_args()

# \section cc Configuration and Initial conditions
#
# We load the Config object from a file
# 
cfg = swarmng.Config.load(args.config)
inputFile = cfg["input"]
snapshotFileName = inputFile + ".snapshot"

#
# Here we check for existence of the snapshot file, if it exists, then we load
# it instead of the initial conditions.
if path.exists(snapshotFileName) :
    fn = snapshotFileName
else :
    fn = inputFile
ext = path.splitext(fn)[1]
if ext == "txt" : ens = swarmng.DefaultEnsemble.load_from_text(fn)
else : ens = swarmng.DefaultEnsemble.load_from_bin(fn)

# \section Resuming the integration
# Using some functional features, we can easily calculate the median of 
# all system times to estimate a starting time for the whole ensemble.
times = sorted(map(lambda s : s.time, ens))
median_of_system_time = times[len(times)/2]
#
# We set the starting time to the median of all system times. 
# We add the duration to the starting time to find a suitable end time.
#
starting_time = median_of_system_time
dt = min(float(cfg["destination_time"]), starting_time + args.duration)
print("Integrating from {0} to {1}".format(starting_time, dt))

#
#  Let the integration begin, it is the same as the @ref TutorialPython "first tutorial"
swarmng.init(cfg)
integ = swarmng.Integrator.create( cfg )
integ.ensemble = ens
integ.destination_time = dt
integ.integrate()

#
# Finally we have to save the snapshot for resuming the integration later
ens.save_to_bin(snapshotFileName)
#
# \section conc Conclusion
#
# In this tutorial, we explored how to save/load snapshot files.
# You can use snapshot files in many integration schemes
# where complex and long simulations are required and
# you may need to end them early sometimes.
#
# Where to go from here:
# * Read documentation for @ref swarmng.DefaultEnsemble for details of load/save methods.
# * Read documentation for @ref swammng.System and @ref swarmng.Body for other attributes on systems and bodies.
# * Follow along to @ref TutorialPythonIntegrationPlot "next tutorial"
