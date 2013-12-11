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


def benchmark_integrator(label,cfg):
    integ = swarmng.Integrator.create(cfg)

    nslist = list(geometric_progression(math.sqrt(2.0),1024,8))
    for nb in range(3,7):
        times = []
        for ns in nslist:
            integ.ensemble = swarmng.generate_ensemble(swarmng.config(nsys=ns,nbod=nb))
            integ.destination_time = 2.0

            t = integ.integrate()
            times.append(t)

        P.plot(nslist,times,label=str(nb));

    P.legend()
    P.savefig("benchmark{0}".format(label))

benchmark_integrator("BPT",swarmng.config(integrator="hermite_bpt_lean",time_step=0.001))
benchmark_integrator("OMP",swarmng.config(integrator="hermite_omp_lean",time_step=0.001))


