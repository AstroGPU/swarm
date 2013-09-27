#!/usr/bin/env python2
# -*- coding: utf8 -*-
## @file logrecord_tutorial.py Tutorial for using the @ref swarmng.logrecord.LogRecord class for accessing records from a BDB log file.

# @page TutorialPythonLogRecord Examining records from the log file.
#
# The Swarm-NG log files consists of records, these records are
# encoded in binary and are accessed through swarmng.logrecord.LogRecord
# class. To get an instance of LogRecord, you would need to run a query
# on a log file and get the records in a for loop, to learn how
# to do that refer to @ref TutorialPythonQueryPlot.
#
# Source code for this tutorial can be found at @ref py/logrecord_tutorial.py

# To run the examples in this tutorial, change the code
# in @ref TutorialPythonQueryPlot to use the function defined 
# here. by default the plot tutorial uses get_semi_major_axis from
# this module.
#
# The easiest way to access the planets in the log record is by
# using bodies_in_keplerian method. This methods converts the
# carteisan coordinates into Keplerian coordinates based on the
# center given. The result is an iterable object that has a 3-tuple
# as elements: 
#  * first: index of the planet, starting from 1 (Note that 0 is the star)
#  * second: the original body attributes in Cartesian coordinates: contains position, velocity and mass.
#  * third: object with orbital elements as properties, refer to \ref swarmng.logrecord.keplerian_for_cartesian for a full list of attributes.
def get_semi_major_axis(l,body):
  for i, b, orbit in l.bodies_in_keplerian(center=l.star()):
    if i == body:
      return orbit.a
