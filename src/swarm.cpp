#include <cuda_runtime_api.h>
#include "swarm.h"
#include <memory>
#include <cmath>
#include <iostream>
#include <valarray>
#include <vector>

// compute the energies of masless bodies in an external potential
// NOTE: This is a temporary function to support the mock "Euler" integrator
void one_body_energy(std::valarray<double> &E, const cpu_ensemble &ens)
{
	E.resize(ens.nbod()*ens.nsys());

	for(int sys = 0; sys != ens.nsys(); sys++)
	{
		for(int bod = 0; bod != ens.nbod(); bod++)
		{
			float m; double x[3], v[3];
			ens.get_body(sys, bod, m, x[0], x[1], x[2], v[0], v[1], v[2]);

			// assuming one-body problem
			E[sys*ens.nbod() + bod] = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) / 2. - 1. / sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
		}
	}
}

// aux class to sort indices by energy error (biggest first)
struct energy_sorter
{
	const std::valarray<double> &dE;
	int nbod;

	energy_sorter(const std::valarray<double> &dE_, int nbod_) : dE(dE_), nbod(nbod_) {}
	bool operator()(const std::pair<int,int>  &a, const std::pair<int,int>  &b) const
	{
		double dEa = fabs(dE[a.first*nbod + a.second]);
		double dEb = fabs(dE[b.first*nbod + b.second]);
		if(isnan(dEa)) { dEa = 0.; }
		if(isnan(dEb)) { dEb = 0.; }
		return dEa > dEb;
	}
};

// just a simple dump to stdout of the first N planets given by indices in idx
void write_output_aux(const std::vector<std::pair<int,int> > &idx, const int N, const cpu_ensemble &ens, const std::valarray<double> &E, const std::valarray<double> &dEoverE)
{
	for(int i = 0; i != N; i++)
	{
		int sys = idx[i].first;
		int bod = idx[i].second;

		float m; double x[3], v[3];
		ens.get_body(sys, bod, m, x[0], x[1], x[2], v[0], v[1], v[2]);

		printf("%5d %5d  T=%f  m=%f  pos=(% 9.5f % 9.5f % 9.5f)  vel=(% 9.5f % 9.5f % 9.5f)  E=% 9.5f dE/E=% 9.5f\n", sys, bod, ens.T(sys), m, x[0], x[1], x[2], v[0], v[1], v[2], E[sys*ens.nbod() + bod], dEoverE[sys*ens.nbod() + bod]);
	}
}

void write_output(cpu_ensemble &ens, std::valarray<double> &E)
{
	std::valarray<double> Enew;
	one_body_energy(Enew, ens);
	std::valarray<double> dEoverE(Enew.size());
	bool has_dE = E.size() == Enew.size();
	if(has_dE)
	{
		dEoverE = Enew / E - 1.;
	}
	else
	{
		E.resize(Enew.size());
		dEoverE = 0.;
	}
	E = Enew;

	std::vector<std::pair<int,int> > idx; idx.reserve(ens.nsys()*ens.nbod());
	for(int i = 0; i != ens.nsys(); i++) for(int j = 0; j != ens.nbod(); j++) idx.push_back(std::make_pair(i, j));

	// print out the first few systems
	write_output_aux(idx, std::min(2, ens.nsys())*ens.nbod(), ens, E, dEoverE);

	if(has_dE)
	{
		// find the planets with worst energy conservation and print them out
		std::sort(idx.begin(), idx.end(), energy_sorter(dEoverE, ens.nbod()));
		std::cerr << "Planets with worst energy conservation (excluding those w. infinite energy):\n";
		write_output_aux(idx, std::min(2, ens.nsys())*ens.nbod(), ens, E, dEoverE);
	}
}

int main()
{
	// load the ensemble
	cpu_ensemble ens;
	load_ensemble("data", ens);
	std::valarray<double> E;
	write_output(ens, E);

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
	write_output(ens, E);

	// both the integrator & the ensembles are automatically deallocated on exit
	// so there's nothing special we have to do here.
	return 0;
}
