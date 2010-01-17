#!/usr/bin/python
''' Similar to peaShooterGen.py, but throw a binary at the planetary system. '''

import math as M
import random as R
#
#
# define basic parameters for ensemble
# keep the primary mass the same for now
# planets can vary in mass, size, location, and number
# this does not consider stability whatsoever
#
# 
nSystems=128
mPrimary=1. # mass of the primary
massMin=.001/330. # 1 earth-mass minimum
minPlanets=2 # keeps these the same
maxPlanets=2
basePlanetMass=0.001 # work on one Jupiter mass.  Variations around that.
minAU=10. # minimum semi-major axis allowed.  If you are running with fixed time steps, be mindful of this setting
maxAU=100.
pert=0.01 # perturbations for other velocities.
HILLS=3. # make sure the planets are separated by this many Hill radii
timeStart=0.
timeEnd=3000000. # time should be given in yr.
incomingR=206265.
maxUnperturbedImpact=1000.
maxBinarySep=1000.
minBinarySep=100.
numObs=1000 # number of observations allowed
ObserveFile="observeTimes.dat"
VelSig=1./30. #(1km/s in code units)
RANDOM_TIMES=0

def getRadius():
	a=0.
	a=maxAU*(2.*R.random()-1.)
	return a+abs(a)/a+minAU

def createObservingFile():
        R.seed()
        obsTimes=[]
        obsTimes.append(timeStart)
        if RANDOM_TIMES:
                for i in xrange(1,numObs):
                        obsTimes.append(R.uniform(timeStart,timeEnd))
                obsTimes.sort()
        else:
                dt=(timeEnd-timeStart)/float(numObs)
                for i in xrange(1,numObs):
                        obsTimes.append(timeStart+i*dt)

        f=open(ObserveFile,"w")
        for i in xrange(numObs):
                f.write(repr(obsTimes[i])+"\n")
        f.close()
        return 0

def getCollision():
	# first use spherical coordinates to find position in sky
	R.seed()
	phi=2.*M.pi*R.random()
	theta=M.pi*(2.*R.random()-1.)
	x=incomingR*M.cos(theta)*M.cos(phi)
	y=incomingR*M.cos(theta)*M.sin(phi)
	z=incomingR*M.sin(theta)
	alpha=100.
	alpha0=maxUnperturbedImpact/incomingR
	while alpha>alpha0:
        	phiPrime=phi+alpha0*(2.*M.sqrt(R.random())-1.)
        	thetaPrime=theta+alpha0*(2.*M.sqrt(R.random())-1.)
		alpha=M.cos(theta)*M.cos(phi)*M.cos(thetaPrime)*M.cos(phiPrime)
		alpha+=M.cos(theta)*M.sin(phi)*M.cos(thetaPrime)*M.sin(phiPrime)
		alpha+=M.sin(thetaPrime)*M.sin(theta)
		alpha=M.acos(alpha)
		print alpha,alpha0
	print "FOUND ONE!!!!!!!!!!!!!!!!!!!"
        VelSigPert=VelSig*(1.+2./3.*(2.*R.random()-1.))
	vx=-VelSigPert*M.cos(thetaPrime)*M.cos(phiPrime)
	vy=-VelSigPert*M.cos(thetaPrime)*M.sin(phiPrime)
	vz=-VelSigPert*M.sin(thetaPrime)
	return x,y,z,vx,vy,vz

def getBinary():
	# first use spherical coordinates to find position in sky
	R.seed()
	phi=2.*M.pi*R.random()
	theta=M.pi*(2.*R.random()-1.)
        dR=R.uniform(minBinarySep,maxBinarySep)
	dx=dR*.5*M.cos(theta)*M.cos(phi)
	dy=dR*.5*M.cos(theta)*M.sin(phi)
	dz=dR*.5*M.sin(theta)
        omega=((mPrimary+mPrimary)/dR**3)**.5
        VelSigPert=.5*dR*omega
	dvx=VelSigPert*M.cos(theta+M.pi*.5)*M.cos(phi)
	dvy=VelSigPert*M.cos(theta+M.pi*.5)*M.sin(phi)
	dvz=VelSigPert*M.sin(theta+M.pi*.5)
	return dx,dy,dz,dvx,dvy,dvz


def main():
	R.seed()
	for i in xrange(nSystems):
		print "Working on system ",i
		nPlanets=R.randint(minPlanets,maxPlanets)
		buffer="data."+repr(i)
		f=open(buffer,"w")

		# write primary first

		f.write(repr(nPlanets+3)+"\n")
        	buffer=repr(mPrimary)+" 0. 0. 0. 0. 0. 0.\n"
		f.write(buffer)
		listx=[]
		for j in xrange(nPlanets):
			mass=max(basePlanetMass*R.random(),massMin)
			x=getRadius()
			listx.append(x)
			OK=0
			if j==0: OK=1
                	while not OK:
				OK=1
				for k in xrange(j):
					minSep=(x**2)**.5*HILLS*(mass/mPrimary/3.)**(1./3.)
					if(abs(x-listx[k])<minSep):
						OK=0
						x=getRadius()
						listx[j]=x
					

			y=0.;z=0.
			vz=x/abs(x)*M.sqrt(mPrimary/abs(x))
			vx=vz*pert*(R.random()*2.-1.)
			vy=vz*pert*(R.random()*2.-1.)
			buffer=repr(mass)+" "+repr(x)+" "+repr(y)+" "+repr(z)+" "+repr(vx)+" "+repr(vy)+" "+repr(vz)+"\n"
			f.write(buffer)
		
		# now find and write collider information
		x,y,z,vx,vy,vz=getCollision()
		dx,dy,dz,dvx,dvy,dvz=getBinary()
                print dvx,vx,dvy,vy,dvz,vz
		buffer=repr(mPrimary)+" "+repr(x+dx)+" "+repr(y+dy)+" "+repr(z+dz)+" "+repr(vx+dvx)+" "+repr(vy+dvy)+" "+repr(vz+dvz)+"\n"
		f.write(buffer)
		buffer=repr(mPrimary)+" "+repr(x-dx)+" "+repr(y-dy)+" "+repr(z-dz)+" "+repr(vx-dvx)+" "+repr(vy-dvy)+" "+repr(vz-dvz)+"\n"
		f.write(buffer)

		f.close()		

	test=createObservingFile()
        if test!=0:
		print "Error when creating observing file."
		return -1
 
	return 0

test=main()
if test==0: print "Program exited sucessfully."
else: "Error in main. Results may be unreliable."
