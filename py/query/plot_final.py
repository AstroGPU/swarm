# -*- coding: utf-8 *-*

#
#
#

from logdb import IndexedLogDB
from log import LogRecord
from query import *
import matplotlib.pyplot as P
from range_type import *

####### Parameters for the plotting algorithm ##########
# Upper bound on number of records to read from log
MAX_RECORDS = 10000
# Name of the file to read the log records from
fileName = "/scratch/hpc/salehqt/kepler37ecc/kep37_log4.db"
# Range of times to accept in the query, universal accepts all times
time_range = Range.interval((0.0,100.0))
# Range of systems to include in the query
system_range = Range.universal() #interval((0,100))
# Types of events to include in the results
event_range = Range.single(1)
# Bodies to record parameters for, only use interval or single
body_range = Range.interval((1,2))
from math import pi
######  Execution Starts here ########################

# Open the input log file
d = IndexedLogDB(fileName)


### Step1: find out how long every system lasted
def get_final_times():
    final_time = {}
#    for i in system_range :
#        for k, l in d.final_conditions(i):
    for i, l in d.final_conditions(system_range):
            #final_time[i] = l.time/(2*pi)
            print i, l.time
    return final_time

def plot_final_times_vs_systemid(final_time):
    P.scatter(final_time.keys(),final_time.values())
    P.axis([system_range.lower(), system_range.upper(), 10, max(final_time.values()) ])
    P.gca().set_yscale('log',basey=10)
    P.savefig("figure_final_time")


final_time = get_final_times()

#for i in final_time:
#    print i, " " , final_time[i] 
#plot_final_times_vs_systemid(final_time)

exit(1)
# Execute a query on the log, the result is a Python iterable
q = d.query(time_range, system_range, event_range)


STARTING_SIZE = 1024
# Set up some lists to hold our numbers for later
body_parameters = {}
for i in body_range:
    body_parameters[i] = zeros(STARTING_SIZE)
times = zeros(STARTING_SIZE)


record_count = 0
# Main loop of processing the records, we use the 
# truncate function to select only the first MAX_RECORDS entries
# from q. The ones that are selected are return in (k, l) tuples
# k is the primary key for the record, l is a logrecord object.
# for most processings we only need l
for k, l in truncate(MAX_RECORDS,q):

    record_count += 1
    recno = record_count - 1
    # double the size of a if we are out of space
    if record_count > times.size :
        new_size = times.size*2
        times = resize(times, new_size)
        for i in body_range:
            body_parameters[i] = resize(body_parameters[i], new_size)

    # Record time of the record for X-axis
    times[recno] = l.time

    # Convert logrecord to keplerian coordinates, the method 
    # bodies_in_keplerian
    # converts the bodies one-by-one. The coordinates are calculated with
    # respect to a center. The possible options are l.star(), l.barycenter()
    # and l.origin()
    # The results can be used in a for loop. This returns triplets as
    #  - Ordinal number of the planet (first planet has number 1)
    #  - Cartesian coordinates of the planet as an object that has properties
    #  position, velocity and mass
    #  - Keplerian coordinates of the planet as an object with properties
    #  a, e, i, O, w, M
    for i, b, orbit in l.bodies_in_keplerian(center=l.star()):
        # lets record major axis of the orbit of the planet in
        # our list
        body_parameters[recno] = orbit.a


times = resize(times, record_count)
for i in body_range:
    body_parameters[i] = resize(body_parameters[i],record_count)


print record_count, " records were processed "
####### Plotting ##############3

# Once all the numbers are gathered in lists we can
# plot them using matplotlib functions

#P.scatter(body_parameters[1],body_parameters[2]);
for i in body_range:
    P.scatter(times,body_parameters[i])
#P.show()
P.savefig("fig2")


