#!/usr/bin/env python2

# @page TutorialPythonEjection Testing for ejection of planets
#
# In this tutorial, we write a test to examine if the stop_on_ejection monitor works
# as expected. It involves generating a test case that causes some planets to eject and
# examining the ensemble after the integration.
#
# The first step is to import the swarmng module. We also use make_test_case from the 
# @ref TutorialPythonEnsemble "Ensemble generation tutorial" to generate our test case.
import swarmng
from ensemble_tutorial import make_test_case
# \subsection candi Configuration and initial conditions
# RMAX is the maximum radius of the planetary system. If the planets is farther than RMAX from
# the origin, we consider it ejected.
RMAX = 10
#
# Configuratino for the integration, we use CPU version of the PEC² Hermite integrator, it is
# by default checks for ejections when the planets get too far from the origin.
cfg = swarmng.config(
        integrator = "hermite_cpu",
        time_step  = 1e-3,
        nogpu      = 1,
        deactivate_on_ejection = 1,
        rmax       = RMAX
        )
# We have to set the destination_time directly on the integrator object
# and it cannot be included in the Config object.
destination_time = 100
#
# We create an ensemble with very close orbits. This ensures that the planets
# will come close at some point and the result would be an ejection.
ens = make_test_case(nsys=20, nbod = 6, spacing_factor=1.01);
#
# \section int Integration
# Same procedure as in @ref TutorialPython. Set-up the integrator parameters and
# call the method @ref swarmng.Integrator.integrate "integrate".
swarmng.init(cfg)
integ = swarmng.Integrator.create( cfg )
integ.ensemble = ens
integ.destination_time = destination_time
integ.integrate()
#
# \section ex Examine
# After the integration finished, we can look into the ensemble to see if
# in fact the integrator has worked as expected.
#
# 
for s in ens:
  for b in s:
    if( b.distance_to_origin() > RMAX ):
      assert(sys.state == -1)
#
#
# 

