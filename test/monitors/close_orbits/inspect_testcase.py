from sys import path
path.append('lib')
from libswarmng_ext import *

from numpy import *
from math import *
from random import uniform
from sys import argv

def norm(l):
    s = 0
    for x in l:
        s += x**2
    return sqrt(s)


ref = DefaultEnsemble.load_from_text(argv[1])
ens = DefaultEnsemble.load_from_text(argv[2])

print "Time, Is enabled, System ID, Body number , Ratio of ejection"
for i in range(0,ens.nsys):
    for j in range(1,ens.nbod):
        ratio  = norm(ens[i][j].pos)/norm(ref[i][j].pos)
        if(ratio > 5):
			print ens[i].time,",", ens[i].is_enabled() ,",",
			print i,",",j,",", ratio
