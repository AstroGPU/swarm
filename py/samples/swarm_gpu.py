from swarmng import *

cfg = Config()

cfg["integrator"] = "hermite_adap"
cfg["time_step_factor"] = "170e-4"
cfg["min_time_step"] = "1e-7"
cfg["max_time_step"] = "1e-2"
cfg["nsys" ] = "16"
cfg["nbod"] = "3"
cfg["time_step"] = "0.001"
cfg["log_writer"] = "null"

ref = generate_ensemble( cfg )
ens = ref.clone()

init(cfg)

integ = GpuIntegrator.create( cfg )

integ.ensemble = ens
integ.destination_time = 1.0


integ.upload_ensemble()

for i in range(100):
	integ.core_integrate()

sync


integ.download_ensemble()

sync


max_deltaE = find_max_energy_conservation_error( ens, ref)
print("Max energy conservation error ", max_deltaE)
