#include <cuda_runtime_api.h>
#include "swarm.h"
#include <memory>

int main()
{
	// these will ultimately be read from a config file
	const int nsys = 1000;
	const int nbod = 2;
	const float dT = 10;	// duration of integration

	// Particular integrator configuration (will also ultimately be loaded from config file)
	config cfg;
	cfg["integrator"] = "gpu_euler";
	cfg["h"] = "0.001";

	// set up a demo ensemble (... config files in the future)
	cpu_ensemble ens(nsys, nbod);
	for(int sys = 0; sys != nbod; sys++)
	{
		for(int bod = 0; bod != nbod; bod++)
		{
			float m = bod == 0 ? 1. : 0.001;
			float vy = bod == 0 ? 0 : 1.;

			ens.set_body(sys, bod,  m,  bod*1., 0., 0.,   0., vy, 0. );
		}
	}

	// set up the integrator
	std::auto_ptr<integrator> integ( create_integrator(cfg["integrator"], cfg) );

	// perform the integration
	gpu_ensemble gpu_ens(ens);
	integ->integrate(ens, dT);

	// download and write dump the results
	ens.copy_from(gpu_ens);
	for(int sys = 0; sys != nbod; sys++)
	{
		for(int bod = 0; bod != nbod; bod++)
		{
			float m; double x[3], v[3];
			ens.get_body(sys, bod, m, x[0], x[1], x[2], v[0], v[1], v[2]);

			printf("%5d %5d  %f  %f %f %f  %f %f %f", sys, bod, m, x[0], x[1], x[2], v[0], v[1], v[2]);
		}
	}

	// both the integrator & the ensembles are automatically deallocated on exit
	// so there's nothing special we have to do here.
	
	return 0;
}
