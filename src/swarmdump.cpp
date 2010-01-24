#include "swarm.h"
#include "swarmio.h"

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

		double dEoverE = (Enew[sys] - Eold[sys]) / Eold[sys];
		printf("%5d %5d  T=%f  m=%f  pos=(% 9.5f % 9.5f % 9.5f)  vel=(% 9.5f % 9.5f % 9.5f)  E=%g  dE/E=%g\n", sys, bod, ens.time(sys), m, x[0], x[1], x[2], v[0], v[1], v[2], Enew[sys], dEoverE);
	}
}

//
// NOTE: This code currently assumes there are two and only two snapshots in the output file
//
int main()
{
	// load the ensemble
	cpu_ensemble ens;

	ens_reader in("output.bin");
	in >> ens;
	unsigned int nprint = std::min(2, ens.nsys());

	// Calculate energy at beginning of integration
	std::valarray<double> Einit(ens.nsys()), Efinal(ens.nsys());
	calc_total_energy(ens, Einit);

	printf("Snapshot #0 (initial conditions):\n");
	for (unsigned int i = 0;i < nprint;++i)
		write_output(ens, i, Einit, Einit);
	printf("\n");

	// find the final snapshot
	int cnt = 0;
	while(in >> ens) { cnt++; }

	// Calculate energy at end of integration
	calc_total_energy(ens, Efinal);

	// store output
	printf("Snapshot #%d (end of simulation)\n", cnt);
	for (unsigned int i = 0;i < nprint;++i)
		write_output(ens, i, Einit, Efinal);
	printf("\n");

	// find systems with worst E conservation
	std::valarray<double> dEoverE = Efinal / Einit - 1.;
	std::vector<int > idx; idx.reserve(ens.nsys());
	for (int i = 0; i != ens.nsys(); i++) idx.push_back(i);
	std::sort(idx.begin(), idx.end(), energy_sorter(dEoverE));
	std::cerr << "Systems with worst energy conservation (excluding those w. infinite energy):\n";
	for (unsigned int i = 0;i < nprint;++i)
		write_output(ens, idx[i], Einit, Efinal);

	return 0;
}
