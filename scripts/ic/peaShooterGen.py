#!/usr/bin/python

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
even_solid_angle=0
nSystems=2048
mPrimary=1. # mass of the primary
massMin=.001/33. # 10 Earth-mass minimum
massMax=.01 # 10 Jupiter-mass maximum
minPlanets=2 # keeps these the same
maxPlanets=2
basePlanetMass=0.001 # work on one Jupiter mass.  Variations around that.
minAU=10. # minimum semi-major axis allowed.  If you are running with fixed time steps, be mindful of this setting
maxAU=30.
pert=0.01 # perturbations for other velocities.
HILLS=2. # make sure the planets are separated by this many Hill radii
timeStart=0.
timeEnd=20000. # time should be given in yr.
incomingR=206265. # AU
maxUnperturbedImpact=1000. # AU
numObs=1000 # number of observations allowed
ObserveFile="observeTimes.dat"
VelSig=1./29.8 #(1km/s in code units)
RANDOM_TIMES=0

def getUniformLog(b0,b1):
        a0=M.log10(b0)
        a1=M.log10(b1)
        a=10.**(R.uniform(a0,a1))
        return a

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
        VelSigPert=VelSig #*(1.+2./3.*(2.*R.random()-1.))
	x=incomingR*M.cos(theta)*M.cos(phi)
	y=incomingR*M.cos(theta)*M.sin(phi)
	z=incomingR*M.sin(theta)
        impact=M.sqrt(maxUnperturbedImpact**2+2.*mPrimary*maxUnperturbedImpact/VelSigPert**2) # assuming reduced mass is mPmP/(mP+mP)
        #impact=maxUnperturbedImpact
	impact*=M.sqrt(R.random())
	alpha0=impact/incomingR
	alpha=1e-3*alpha0 # take advantage of large separation
	ppf=(2.*(R.random())-1.)*alpha
	tpf=(2.*(R.random())-1.)*alpha
	while alpha<alpha0:
			
       		phiPrime=phi+ppf
               	thetaPrime=theta+tpf
		ppf*=1.01;tpf*=1.01

		alpha=M.cos(theta)*M.cos(phi)*M.cos(thetaPrime)*M.cos(phiPrime)
		alpha+=M.cos(theta)*M.sin(phi)*M.cos(thetaPrime)*M.sin(phiPrime)
		alpha+=M.sin(thetaPrime)*M.sin(theta)
		alpha=M.acos(alpha)
	print alpha,alpha0,impact,alpha*incomingR
	#print "FOUND ONE!!!!!!!!!!!!!!!!!!!"
	vx=-VelSigPert*M.cos(thetaPrime)*M.cos(phiPrime)
	vy=-VelSigPert*M.cos(thetaPrime)*M.sin(phiPrime)
	vz=-VelSigPert*M.sin(thetaPrime)
	return x,y,z,vx,vy,vz

def main():
	R.seed()
	for i in xrange(nSystems):
		#print "Working on system ",i
		nPlanets=R.randint(minPlanets,maxPlanets)
		buffer="data."+repr(i)
		f=open(buffer,"w")

		# write primary first

		f.write(repr(nPlanets+2)+"\n")
        	#buffer=repr(mPrimary)+" 0. 0. 0. 0. 0. 0.\n"
		#f.write(buffer)
                buffer=""
		listx=[]
		pvx=0.;pvy=0.;pvz=0.
		for j in xrange(nPlanets):
			mass=getUniformLog(massMin,massMax)
			x=getUniformLog(minAU,maxAU)
			if j==0: x=10.
			else: x=30.
			listx.append(x)
			OK=0
			if j==0: OK=1
                	while not OK:
				OK=1
				for k in xrange(j):
					minSep=abs(x)*HILLS*(mass/mPrimary/3.)**(1./3.)
					if( abs(abs(x)-abs(listx[k]))<minSep):
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
		
		# now find and write collider information
		x,y,z,vx,vy,vz=getCollision()
		buffer=repr(mPrimary)+" "+repr(x)+" "+repr(y)+" "+repr(z)+" "+repr(vx)+" "+repr(vy)+" "+repr(vz)+"\n"
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
