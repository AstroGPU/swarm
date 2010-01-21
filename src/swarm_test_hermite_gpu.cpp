#include "swarm.h"
#include <memory>
#include <cmath>
#include <iostream>
#include <valarray>
#include <sstream>
#include <fstream>
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

// Convert a variable of arbitrary type to a string.
// NOTE: heavy (unoptimized) function, use sparingly
template<typename T>
std::string str(const T& var)
{
	std::ostringstream ss;
	ss << var;
	return ss.str();
}

// trim whitespaces from the beginning and the end of a string
void trim(std::string& str)
{
	std::string::size_type pos = str.find_last_not_of(" \t");
	if (pos != std::string::npos)
	{
		str.erase(pos + 1);
		pos = str.find_first_not_of(" \t");
		if (pos != std::string::npos) str.erase(0, pos);
	}
	else str.erase(str.begin(), str.end());
}

// load a configuration file
void load_config(config &cfg, const std::string &fn)
{
	std::ifstream in(fn.c_str());
	if(!in) ERROR("Cannot open configuration file '" + fn + "'.");

	std::string line;
	int iline = 0;
	while(std::getline(in, line))
	{
		iline++;
		trim(line);
		if(line.empty()) { continue; }
		if(line[0] == '#') { continue; }

		size_t eqpos = line.find(' ');
		if(eqpos == std::string::npos) ERROR("Error on line " + str(line) + ": '=' sign expected.");

		std::string key = line.substr(0, eqpos), val = line.substr(eqpos+2);
		trim(key); trim(val);

		cfg[key] = val;
	}
}

// get a configuration value for 'key', throwing an error if it doesn't exist
// NOTE: heavy (unoptimized) function, use sparingly
template<typename T>
void get_config(T &val, const config &cfg, const std::string &key)
{
	if(!cfg.count(key)) { ERROR("Configuration key '" + key + "' missing."); }
	std::istringstream ss(cfg.at(key));
	ss >> val;
}

#define SWATCH_STOP(s)  { cudaThreadSynchronize(); (s).stop(); }
#define SWATCH_START(s) { (s).start(); }

int main()
{
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

	// performance stopwatches
	stopwatch swatch_kernel, swatch_mem, swatch_temps, swatch_all;

	// set up the integrator
	SWATCH_START(swatch_all);
	config cfg;
	load_config(cfg, "integrator.cfg");
	std::auto_ptr<integrator> integ(integrator::create(cfg));
	std::string runon = cfg.count("runon") ? cfg["runon"] : "gpu";
	bool ongpu;
	     if(runon == "gpu") { ongpu = true; }
	else if(runon == "cpu") { ongpu = false; }
	else { ERROR("The 'runon' configuration file parameter must be one of 'gpu' or 'cpu'"); }
	std::cerr << "Integrator: " << cfg["integrator"] << ", executing on the " << (ongpu ? "GPU" : "CPU") << "\n";
	// duration of integration
	double dT;
	get_config(dT, cfg, "dT");

	// perform the integration
	if(ongpu)
	{
		SWATCH_START(swatch_mem);
		gpu_ensemble gpu_ens(ens);				// upload to GPU
		SWATCH_STOP(swatch_mem);

		SWATCH_START(swatch_temps);
		integ->integrate(gpu_ens, 0.);				// initialize internal data structures
		SWATCH_STOP(swatch_temps);

		SWATCH_START(swatch_kernel);
		integ->integrate(gpu_ens, dT);				// integrate
		SWATCH_STOP(swatch_kernel);

		SWATCH_START(swatch_mem);
		ens.copy_from(gpu_ens);					// download to host
		SWATCH_STOP(swatch_mem);
	}
	else
	{
		SWATCH_START(swatch_temps);
		integ->integrate(ens, 0.);				// initialize internal data structures
		SWATCH_STOP(swatch_temps);

		SWATCH_START(swatch_kernel);
		integ->integrate(ens, dT);				// integrate
		SWATCH_STOP(swatch_kernel);
	}
	SWATCH_STOP(swatch_all);

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
	std::cerr << "Time per system (integration)   : " << us_per_sys_kernel << " us.\n";
	std::cerr << "Time per system (setup+integr.) : " << us_per_sys_all << " us.\n";
	std::cerr << "GPU/CPU memcpy time             : " << swatch_mem.getTime()*1000 << " ms.\n";
	std::cerr << "Internal state initialization   : " << swatch_temps.getTime()*1000 << " ms.\n";

	// both the integrator & the ensembles are automatically deallocated on exit
	// so there's nothing special we have to do here.
	return 0;
}
