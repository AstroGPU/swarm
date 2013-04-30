from sys import path
path.append('lib')
from libswarmng_ext import *

from numpy import *
from math import *
from random import uniform

cfg = Config()

cfg["integrator"] = "hermite_adap"
cfg["time_step_factor"] = "170e-4"
cfg["min_time_step"] = "1e-7"
cfg["max_time_step"] = "1e-2"


def make_test_case(nsys = 16, nbod = 3 , spacing_factor = 1.4, planet_mass = 0.001, ejection_factor = 1):
    d = DefaultEnsemble.create(nbod,nsys)
    for i in range(0,d.nsys):
        s = d[i]
        s.id = i
        s.set_active()
        s.time = 0
        s[0].pos = [ 0, 0, 0 ]
        s[0].vel = [ 0, 0, 0 ]
        s[0].mass = 1
        for j in range(1,d.nbod):
            r = spacing_factor ** (j-1)
            v = sqrt(1/r) * ejection_factor
            phi = uniform(0,2*pi)
            s[j].pos = [  r*cos(phi), r*sin(phi), 0 ]
            s[j].vel = [ -v*sin(phi), v*cos(phi), 0 ]
            s[j].mass = planet_mass 
    return d;


def norm(l):
    s = 0
    for x in l:
        s += x**2
    return sqrt(s)


# Initializations
init(cfg)
integ = Integrator.create( cfg )


# Integrating
ref = make_test_case(nsys=20, nbod = 6, spacing_factor=1.01);

ref.save_to_text("hh.txt")

ens = ref.clone()
integ.ensemble = ens
integ.destination_time = 100.0
integ.integrate()

for i in range(0,ens.nsys):
    for j in range(0,ens.nbod):
        ratio  = norm(ens[i][j].pos)/norm(ref[i][j].pos)
        if(ratio > 5):
            print i,j, ratio

