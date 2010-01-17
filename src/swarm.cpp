#include <cuda_runtime_api.h>
#include "swarm.h"
#include <memory>
#include <cmath>

void load_ensemble(cpu_ensemble &ens)
{
	// set up a demo ensemble (TODO: load from config files in the future)
	const int nsys = 1000;
	const int nbod = 2;
	ens.reset(nsys, nbod);

	// generate nbod bodies with random positions and velocities,
	const double twopi = 3.141592653589793;
	srand(42);
	for(int sys = 0; sys != nsys; sys++)
	{
		for(int bod = 0; bod != nbod; bod++)
		{
			float m = 0.001;
			float x  =  2.f * ( ((float)random()) / RAND_MAX - 0.5 );
			float y  =  2.f * ( ((float)random()) / RAND_MAX - 0.5 );
			float z  =  2.f * ( ((float)random()) / RAND_MAX - 0.5 );
			float vx = 0.5f * ( ((float)random()) / RAND_MAX - 0.5 );
			float vy = 0.5f * ( ((float)random()) / RAND_MAX - 0.5 );
			float vz = 0.5f * ( ((float)random()) / RAND_MAX - 0.5 );

// 			float a = 10.f * ( ((float)random()) / RAND_MAX ); // a = [0..10)
// 			float e = 1.f * ( ((float)random()) / RAND_MAX ); // e = [0..1)
// 			float i = twopi * 0.5 * ( ((float)random()) / RAND_MAX );
// 			float o = twopi * ( ((float)random()) / RAND_MAX );
// 			float O = twopi * ( ((float)random()) / RAND_MAX );
// 			float M = twopi * ( ((float)random()) / RAND_MAX );

			ens.set_body(sys, bod, m, x, y, z, vx, vy, vz);
		}
	}
}

void write_output(cpu_ensemble &ens)
{
	// just a simple dump to stdout
	for(int sys = 0; sys != std::min(5, ens.nsys()); sys++)
	{
		for(int bod = 0; bod != ens.nbod(); bod++)
		{
			float m; double x[3], v[3];
			ens.get_body(sys, bod, m, x[0], x[1], x[2], v[0], v[1], v[2]);

			// assuming one-body problem
			float E = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) / 2. - 1. / sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

			printf("%5d %5d  T=%f  m=%f  pos=(%f %f %f)  vel=(%f %f %f)  E=%f\n", sys, bod, ens.T(sys), m, x[0], x[1], x[2], v[0], v[1], v[2], E);
		}
	}
}

int main()
{
	// load the ensemble
	cpu_ensemble ens;
	load_ensemble(ens);
	write_output(ens);

	// set up the integrator and integrator config (TODO: load from config file)
	config cfg;
	cfg["integrator"] = "gpu_euler";
	//cfg["integrator"] = "gpu_mvs";
	//cfg["integrator"] = "gpu_does_not_exist";
	cfg["h"] = "0.001";
	std::auto_ptr<integrator> integ( integrator::create(cfg) );

	// perform the integration
	const float dT = 1;	// duration of integration (TODO: read from config file)
	gpu_ensemble gpu_ens(ens);
	integ->integrate(gpu_ens, dT);
	ens.copy_from(gpu_ens);

	// store output
	write_output(ens);

	// both the integrator & the ensembles are automatically deallocated on exit
	// so there's nothing special we have to do here.
	return 0;
}
