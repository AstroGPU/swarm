from math import *
from random import uniform
import swarmng
import unittest

##  Testing the monitor stop_on_ejection, some planets are put in orbits
##  very close to one another (1.01 separation) and they are expected to
##  eject out of the system. Since not all of the eject, in the
##  verification stage we check that if the body has ejected, the system
##  should have been deactivated

def make_test_case(nsys = 16, nbod = 3 , spacing_factor = 1.4, planet_mass = 0.001, ejection_factor = 1):
	d = swarmng.DefaultEnsemble.create(nbod,nsys)
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

class EjectionTest(unittest.TestCase):
	RMAX = 10
	def runTest(self):
		cfg = swarmng.mkConfig({
		  'integrator' : "hermite_cpu",
		  'time_step' : 170e-4,
		  "min_time_step" : 1e-7,
		  'max_time_step': 1e-2 ,
		  'nogpu' : 1,
		  'deactivate_on_ejection': '1',
		  'rmax' : self.RMAX,
		  })

		# Initializations
		swarmng.init(cfg)
		integ = swarmng.Integrator.create( cfg )

		# Integrating
		ref = make_test_case(nsys=20, nbod = 6, spacing_factor=1.01);

		ens = ref.clone()
		integ.ensemble = ens
		integ.destination_time = 100.0
		integ.integrate()

		for sys in ens:
			for b in sys:
				if( b.distance_to_origin() > self.RMAX ):
					self.assertEqual(sys.state, -1, "a system with an ejected body should be disabled")



