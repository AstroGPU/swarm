#!/usr/bin/python
#
#    "throwBinaryGen.py" is a python script that creates initial conditions and the observation
#    file for use in swarm_scatter_demo. Like, swarm_scatter_demo.py, except it throws binaries
#    instead of single perturbers.
#
#    "swarm_scatter_demo" is a program that uses the Swarm-NG tools for modeling an ensemble of
#    small N systems using the hermite_adap_gpu integrator.
#    Copyright (C) 2010  Swarm-NG Development Group
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Contact aaron.boley at the domain gmail.com if you questions regarding
#    this software.
#
#
import math as M
import random as R
import scipy.special as SS
#
#
# define basic parameters for ensemble
# keep the primary mass the same for now, and binary masses are assumedto be equal
# planets can vary in mass, size, location, and number
# this does not consider stability, except for the Hill radius check
# the outer planet is always set to 128 AU, but that can be changed where indicated
# perturbers are drawn from a Maxwellian velocity distribution
#
# 
nSystems=2048 # systems in ensemble
mPrimary=1. # mass of the primary
mBinary=mPrimary/2. # mass of binary
nOther=3 # primary plus perturbers. E.g., 3 = primary plus one binary
massMin=.001/32. # 10 Earth-mass minimum
massMax=.001 # 1 Jupiter-mass maximum
minPlanets=4 # keeps these the same for now.
maxPlanets=4
minAU=2. # minimum semi-major axis allowed.  If you are running with fixed time steps, be mindful of this setting.
maxAU=100. # outer planet
minBinarySep=100. # binary separations
maxBinarySep=100.
pert=0.01 # perturbations for other velocities.
HILLS=3. # make sure the planets are separated by HILLS many Hill radii
timeStart=0.
timeEnd=2e5 # time should be given in yr.
incomingR=20626.5 # AU
maxUnperturbedImpact=1000. # AU
numObs=2000 # number of observations allowed
ObserveFile="observeTimes.dat"
VelSig=1./29.8 #(1km/s in code units) This is mean of Maxwellian velocity dispersion
RANDOM_TIMES=0
thisSeed=314159

def getUniformLog(b0,b1):
        a0=M.log10(b0)
        a1=M.log10(b1)
        a=10.**(R.uniform(a0,a1))
        return a

def Maxwell(v):
	p=R.random()
        s=0.
        SqrtPi=M.sqrt(M.pi)
        xinit=.001
        x=xinit
        while True:
                s=M.pi/v*(SqrtPi*v*SS.erf(x/v)-2*x*M.exp(-(x/v)**2))/M.pi**1.5
                err=(s-p)/p
                if M.fabs(err)<1e-6:break
                if(err<0.):x=x-x*err*.1
                else:x=x-x*err*.1
                if x<0.:
                        xinit/=2.
                        x=xinit
        return x

def createObservingFile():
        obsTimes=[]
        obsTimes.append(timeStart)
        if RANDOM_TIMES:
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

def getCollision():
	# first use spherical coordinates to find position in sky as observed by planetary system
	phi=2.*M.pi*R.random()
	theta=M.pi*(2.*R.random()-1.)
        #VelSigPert=Maxwell(VelSig) 
	VelSigPert=VelSig
	x=incomingR*M.cos(theta)*M.cos(phi)
	y=incomingR*M.cos(theta)*M.sin(phi)
	z=incomingR*M.sin(theta)
        impact=M.sqrt(maxUnperturbedImpact**2+2.*mBinary*maxUnperturbedImpact/VelSigPert**2) # assuming reduced mass is mPmP/(mP+mP)
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
	vx=-VelSigPert*M.cos(thetaPrime)*M.cos(phiPrime)
	vy=-VelSigPert*M.cos(thetaPrime)*M.sin(phiPrime)
	vz=-VelSigPert*M.sin(thetaPrime)
	return x,y,z,vx,vy,vz

def getBinary():
	phi=2.*M.pi*R.random()
	theta=M.pi*(2.*R.random()-1.)
        dR=R.uniform(minBinarySep,maxBinarySep)
	dx=dR*.5*M.cos(theta)*M.cos(phi)
	dy=dR*.5*M.cos(theta)*M.sin(phi)
	dz=dR*.5*M.sin(theta)
        omega=((mBinary+mBinary)/dR**3)**.5
        VelSigPert=.5*dR*omega
	dvx=VelSigPert*M.cos(theta+M.pi*.5)*M.cos(phi)
	dvy=VelSigPert*M.cos(theta+M.pi*.5)*M.sin(phi)
	dvz=VelSigPert*M.sin(theta+M.pi*.5)
	return dx,dy,dz,dvx,dvy,dvz

def main():
	R.seed(thisSeed)
	for i in xrange(nSystems):
		nPlanets=R.randint(minPlanets,maxPlanets)
		print " Working on system ", i
		buffer="data."+repr(i)
		f=open(buffer,"w")

		f.write(repr(nPlanets+nOther)+"\n")
                buffer=""
		listx=[]
		pvx=0.;pvy=0.;pvz=0.
		for j in xrange(nPlanets):
			mass=getUniformLog(massMin,massMax)
			x=getUniformLog(minAU,maxAU)
			if j==0: x=maxAU
			if j==1: x=3.
			if j==2: x=10. # set these if you want to constrain planet locations
			if j==3: x=30.
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
		
		# now find and write perturber information
		for j in xrange(nOther-2):
			x,y,z,vx,vy,vz=getCollision()
			dx,dy,dz,dvx,dvy,dvz=getBinary()
			buffer=repr(mBinary)+" "+repr(x+dx)+" "+repr(y+dy)+" "+repr(z+dz)+" "+repr(vx+dvx)+" "+repr(vy+dvy)+" "+repr(vz+dvz)+"\n"
			f.write(buffer)
			buffer=repr(mBinary)+" "+repr(x-dx)+" "+repr(y-dy)+" "+repr(z-dz)+" "+repr(vx-dvx)+" "+repr(vy-dvy)+" "+repr(vz-dvz)+"\n"
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
