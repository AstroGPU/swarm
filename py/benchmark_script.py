#!/usr/bin/env python2
# -*- coding: utf8 -*-

import swarmng
import time
import matplotlib.pyplot as P

import math

def geometric_progression(base,init,count):
    q = base
    p = init
    for i in range(0,count):
        yield p
        p = q * p

swarmng.init(swarmng.config(verbose=1))

integrators = [ "hermite_bpt_lean", "hermite_omp_lean" ]

def benchmark_integrator(nb,cfg):
    """For one specific number of bodies compare different
    integrators, show a plot and write the results to
    a file"""


    nslist = list(geometric_progression(math.sqrt(2.0),1024,15))
    result = {"nsys" : nslist, "nbod": nb }
    for I in integrators:

        cfg["integrator"] = I
        integ = swarmng.Integrator.create(cfg)

        times = []; times_map = {}
        for ns in nslist:
            print("Integrating {0} {1} {2}".format(nb,I,ns))
            integ.ensemble = swarmng.generate_ensemble(swarmng.config(nsys=ns,nbod=nb))
            integ.destination_time = 2.0

            t = integ.integrate()
            times.append(t)

        result[I] = times
        #P.plot(nslist,times,label=I);

    #P.legend()
    #P.show()
    return result

r = benchmark_integrator(3,swarmng.config(time_step=0.001))

import json
print(json.dumps(r, sort_keys = True, indent=2))


