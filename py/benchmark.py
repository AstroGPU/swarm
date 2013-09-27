#!/usr/bin/env python2
# -*- coding: utf8 -*-

import swarmng
import time

for c in range(0,2):
    swarmng.init(swarmng.config(CUDA_DEVICE=c,verbose=1))

    integ = swarmng.Integrator.create(swarmng.config(integrator="hermite",time_step=0.001))

    times = []

    for nb in range(3,7):
        integ.ensemble = swarmng.generate_ensemble(swarmng.config(nsys=8000,nbod=nb))
        integ.destination_time = 10.0

        start = time.clock()
        integ.integrate()
        times.append(time.clock() - start)

    print(times)




