# Header
# header 


# @page TutorialPython Beginner Python Tutorial
import sys
import os

BASEDIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append('lib')

# raoecuh
#a.,rcuhc,r.uh
from numpy import *
from libswarmng_ext import *
# ar,c.uh.,cu
# ourcha

cfg = Config()

cfg["integrator"] = "hermite_cpu"
cfg["time_step_factor"] = "170e-4"
cfg["min_time_step"] = "1e-7"
cfg["max_time_step"] = "1e-2"
cfg["nsys" ] = "16"
cfg["nbod"] = "3"
cfg["time_step"] = "0.001"
cfg["log_writer"] = "null"
cfg["nogpu"] = "1"


def integrate(ens) :
	init(cfg)
	integ = Integrator.create( cfg )
	integ.ensemble = ens
	integ.destination_time = 1.0
	sync
	integ.integrate()
	sync


ref = DefaultEnsemble.load_from_text( os.path.join(BASEDIR , "../test/integrators/test.3.in.txt"))

for s in range(ref.nsys):
	for b in range(ref.nbod):
		print("pos: %s\tvel:%s" % (ref[s][b].pos, ref[s][b].vel))

ref[0][0].pos = [ 1.0, 2.0 , 5.0 ]

ens = ref.clone()
integrate(ens)

max_deltaE = find_max_energy_conservation_error( ens, ref)
print("Max energy conservation error %g" % max_deltaE)

