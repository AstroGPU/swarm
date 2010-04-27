#!/usr/bin/python
''' Simple python script for generating N systems with M planets '''

import math as M
import random as R
#
#
# define basic parameters for ensemble
# keep the primary mass the same for now
# planets can vary in mass, size, location, and number
# this does not consider stability execpt for simple
# Hill radius criterion
#
# Example parameters 
#
nSystems=1024
mPrimary=1. # mass of the primary
massMin=.001/320. # 1 earth-mass minimum
massMax=0.01 # 10 Jupiter mass max
minPlanets=2 # keeps these the same
maxPlanets=2
minAU=1.0 # minimum semi-major axis allowed.  If you are running with fixed time steps, be mindful of this setting
maxAU=10.
pert=0.01 # perturbations for other velocities
HILLS=3. # make sure the planets are separated by this many Hill radii
timeStart=0.
timeEnd=100. # time should be given in yr
numObs=100 # number of system observations. File is unnecessary for most demos. Used as the list of times file in swarm_scatter_demo
ObserveFile="observeTimes.dat"
RANDOM_TIMES=0
thisSeed = 314159 # Seed for random generator

def getUniformLog(b0,b1):
        a0=M.log10(b0)
        a1=M.log10(b1)
        a=10.**(R.uniform(a0,a1))
        return a

def createObservingFile():
	obsTimes=[]
        obsTimes.append(timeStart) 
        if RANDOM_TIMES:
		R.seed(thisSeed)
		for i in xrange(1,numObs):
			obsTimes.append(R.uniform(timeStart,timeEnd))
		obsTimes.sort()
        else:
		dt=(timeEnd-timeStart)/float(numObs)
		for i in xrange(1,numObs+1):
			obsTimes.append(timeStart+i*dt)

        f=open(ObserveFile,"w")
	for i in xrange(numObs+1):
		f.write(repr(obsTimes[i])+"\n")
	f.close()
	return 0

def main():
	R.seed(thisSeed)
	for i in xrange(nSystems):
		# print "Working on system ",i
		nPlanets=R.randint(minPlanets,maxPlanets)
		buffer="data."+repr(i)
		f=open(buffer,"w")
		buffer=""

		# write primary first

		f.write(repr(nPlanets+1)+"\n")
		listx=[]
		pvx=0.;pvy=0.;pvz=0.
		for j in xrange(nPlanets):
                        mass=getUniformLog(massMin,massMax)
                        x=getUniformLog(minAU,maxAU)
			listx.append(x)
			OK=0
			if j==0: OK=1
                	while not OK:
				OK=1
				for k in xrange(j):
					minSep=(x**2)**.5*HILLS*(mass/mPrimary/3.)**(1./3.)
					if(abs(abs(x)-abs(listx[k]))<minSep):
						OK=0
                                                x=getUniformLog(minAU,maxAU)
						listx[j]=x
					
			y=0.;z=0.
			vy=x/abs(x)*M.sqrt(mPrimary/abs(x))
			vx=vy*pert*(R.random()*2.-1.)
			vz=vy*pert*(R.random()*2.-1.)
                        pvy-=vy*mass
			pvx-=vx*mass
			pvz-=vz*mass
			buffer+=repr(mass)+" "+repr(x)+" "+repr(y)+" "+repr(z)+" "+repr(vx)+" "+repr(vy)+" "+repr(vz)+"\n"
		buffer0=repr(mPrimary)+" "+repr(0.)+" "+repr(0.)+" "+repr(0.)+" "+repr(pvx/mPrimary)+" "+repr(pvy/mPrimary)+" "+repr(pvz/mPrimary)+"\n"
		f.write(buffer0)
		f.write(buffer)
		f.close()		

	print "Generated",i+1,"systems."
	test=createObservingFile()
        if test!=0:
		print "Error when creating observing file."
		return -1
 
	return 0

test=main()
if test!=0: print "Error in main. Results may be unreliable."
