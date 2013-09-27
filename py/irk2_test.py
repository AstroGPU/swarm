import swarmng, argparse
import os

swarmng.init(swarmng.config(nogpu=0))

ensemble_file = os.path.join(swarmng.SWARMDIR,'test/integrators/test.3.in.txt')

ref = swarmng.DefaultEnsemble.load_from_text(ensemble_file)

ens = ref.clone()

integrator_cfg = swarmng.config(
	integrator = 'irk2',
	time_step = .001,
	nbod = 3
)

integ = swarmng.Integrator.create(integrator_cfg)

integ.ensemble = ens
integ.destination_time = 100

integ.integrate()

ens.save_to_text("hermite_cpu_snapshot.txt")

max_deltaE = swarmng.find_max_energy_conservation_error(ens,ref)
print("Max energy conservation error %g" % max_deltaE)
