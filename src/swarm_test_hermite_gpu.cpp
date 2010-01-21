#include "swarm.h"
#include <memory>
#include <cmath>
#include <iostream>
#include <valarray>
#include <vector>

// aux class to sort indices by energy error (biggest first)
struct energy_sorter
{
	const std::valarray<double> &dE;

	energy_sorter(const std::valarray<double> &dE_) : dE(dE_) {};
	bool operator()(const int  &a, const int  &b) const
	{
		double dEa = fabs(dE[a]);
		double dEb = fabs(dE[b]);
		if (isnan(dEa)) { dEa = 0.; }
		if (isnan(dEb)) { dEb = 0.; }
		return dEa > dEb;
	}
};

// just a simple dump to stdout of one system
void write_output(const cpu_ensemble &ens, const int sys, std::valarray<double>  &Eold, std::valarray<double> &Enew)
{
	for (int bod = 0; bod != ens.nbod(); bod++)
	{
		float m; double x[3], v[3];
		ens.get_body(sys, bod, m, x[0], x[1], x[2], v[0], v[1], v[2]);

		printf("%5d %5d  T=%f  m=%f  pos=(% 9.5f % 9.5f % 9.5f)  vel=(% 9.5f % 9.5f % 9.5f)  E=%g", sys, bod, ens.time(sys), m, x[0], x[1], x[2], v[0], v[1], v[2], Enew[sys]);
		if (Enew[sys] != 0)
			printf("  dE/E=%g\n", (Enew[sys] - Eold[sys]) / Eold[sys]);
		else
			printf("\n");
	}
}

int main()
{
	// perform the integration
	stopwatch swatch_kernel, swatch_all;
	swatch_all.start();

	// load the ensemble
	cpu_ensemble ens;
	load_ensemble("data", ens);
	unsigned int nprint = std::min(2, ens.nsys());


	// Calculate energy at beginning of integration
	std::valarray<double> Einit(ens.nsys()), Efinal(ens.nsys());
	calc_total_energy(ens, Einit);

	printf("Initial conditions...\n");
	for (unsigned int i = 0;i < nprint;++i)
		write_output(ens, i, Einit, Efinal);

	// set up the integrator and integrator config (TODO: load from config file)
	config cfg;
	cfg["integrator"] = "gpu_hermite";
	//cfg["integrator"] = "cpu_hermite";
	//cfg["integrator"] = "gpu_mvs";
	//cfg["integrator"] = "gpu_does_not_exist";
	cfg["h"] = "0.001";
	cfg["precision"] = "1";
	std::auto_ptr<integrator> integ(integrator::create(cfg));

	// perform the integration
	const float dT = 1;	// duration of integration (TODO: read from config file)

	gpu_ensemble gpu_ens(ens);				// upload to GPU
	swatch_kernel.start();
	integ->integrate(gpu_ens, dT);				// integrate
	//integ->integrate(ens, dT);
	cudaThreadSynchronize(); swatch_kernel.stop();
	ens.copy_from(gpu_ens);					// download to host
	cudaThreadSynchronize(); swatch_all.stop();

	// Calculate energy at end of integration
	calc_total_energy(ens, Efinal);

	// store output
	printf("Final conditions...\n");
	for (unsigned int i = 0;i < nprint;++i)
		write_output(ens, i, Einit, Efinal);

	// find systems with worst E conservation
	std::valarray<double> dEoverE = Efinal / Einit - 1.;
	std::vector<int > idx; idx.reserve(ens.nsys());
	for (int i = 0; i != ens.nsys(); i++) idx.push_back(i);
	std::sort(idx.begin(), idx.end(), energy_sorter(dEoverE));
	std::cerr << "Systems with worst energy conservation (excluding those w. infinite energy):\n";
	for (unsigned int i = 0;i < nprint;++i)
		write_output(ens, idx[i], Einit, Efinal);

	// print out timings
	double us_per_sys_all = (swatch_all.getTime() / ens.nsys()) * 1000000;
	double us_per_sys_kernel = (swatch_kernel.getTime() / ens.nsys()) * 1000000;
	std::cerr << "Time per system (integration): " << us_per_sys_kernel << " us.\n";
	std::cerr << "Time per system (setup+integr.): " << us_per_sys_all << " us.\n";

	// both the integrator & the ensembles are automatically deallocated on exit
	// so there's nothing special we have to do here.
	return 0;
}
