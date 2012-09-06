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


ens = DefaultEnsemble.load_from_text(argv[1])

print "Time, Body number 1, Body number 2, Distance between bodies"
for i in range(0,ens.nsys):
	if ens[i].state == -1:
		print ens[i].time,
		for j in range(1,ens.nbod):
			for k in range(1,j):
				a = ens[i][j].pos
				b = ens[i][k].pos
				dist = sqrt((a[0]-b[0])**2 +(a[1]-b[1])**2+(a[2]-b[2])**2)
				if dist < .01:
					print j,k, dist,
		print ""
