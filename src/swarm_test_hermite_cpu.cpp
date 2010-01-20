#include "swarm.h"
#include <memory>
#include <cmath>
#include <iostream>
#include <valarray>
#include <vector>

// compute the energies of each system
// NOTE: This is different from the one used for Euler integrator
double calc_system_energy(const cpu_ensemble &ens, const int sys)
{
  double E = 0.;
  
/*  for(int sys = 0; sys != ens.nsys(); sys++)
    {*/
//      E = 0.;
      for(int bod1 = 0; bod1 != ens.nbod(); bod1++)
	{
	  float m1; double x1[3], v1[3];
	  ens.get_body(sys, bod1, m1, x1[0], x1[1], x1[2], v1[0], v1[1], v1[2]);
	  E += 0.5*m1*(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
	  
	  for(int bod2 = 0; bod2 < bod1; bod2++)
	    {
	      float m2; double x2[3], v2[3];
	      ens.get_body(sys, bod2, m2, x2[0], x2[1], x2[2], v2[0], v2[1], v2[2]);
	      double dist = sqrt((x2[0]-x1[0])*(x2[0]-x1[0])+(x2[1]-x1[1])*(x2[1]-x1[1])+(x2[2]-x1[2])*(x2[2]-x1[2]));
	      
	      E -= m1*m2/dist;
	    }
	}
	return E;
/*    }*/
}

// just a simple dump to stdout of one system
void write_output(const cpu_ensemble &ens, const int sys, double& Eold)
{
  double Enew = calc_system_energy(ens,sys);
  for(int bod = 0; bod != ens.nbod(); bod++)
    {
      float m; double x[3], v[3];
      ens.get_body(sys, bod, m, x[0], x[1], x[2], v[0], v[1], v[2]);
      
      printf("%5d %5d  T=%f  m=%f  pos=(% 9.5f % 9.5f % 9.5f)  vel=(% 9.5f % 9.5f % 9.5f)  E=%g", sys, bod, ens.time(sys), m, x[0], x[1], x[2], v[0], v[1], v[2],Enew);
      if(Eold!=0)
	printf("  dE/E=%g\n",(Enew-Eold)/Eold);
      else
	printf("\n");
    }
  Eold = Enew;
}

int main()
{
	// load the ensemble
	cpu_ensemble ens;
	load_ensemble("data", ens);
	std::vector<double> E(ens.nsys(),0.);

	printf("Initial conditions...\n");
	write_output(ens, 0, E[0]);
	write_output(ens, 1, E[1]);
	write_output(ens, 2, E[2]);

	// set up the integrator and integrator config (TODO: load from config file)
	config cfg;
	cfg["integrator"] = "cpu_hermite";
	//cfg["integrator"] = "gpu_mvs";
	//cfg["integrator"] = "gpu_does_not_exist";
	cfg["h"] = "0.001";
	std::auto_ptr<integrator> integ( integrator::create(cfg) );

	// perform the integration
	const float dT = 1;	// duration of integration (TODO: read from config file)

	/* Replace GPU call here with CPU call below comment
	   gpu_ensemble gpu_ens(ens);
	   integ->integrate(gpu_ens, dT);
	   ens.copy_from(gpu_ens);
	*/


//	integ->allocate_workspace_for_ensemble(ens);
	integ->integrate(ens, dT);

	// store output
	printf("Final conditions...\n");
	write_output(ens, 0, E[0]);
	write_output(ens, 1, E[1]);
	write_output(ens, 2, E[2]);

	// both the integrator & the ensembles are automatically deallocated on exit
	// so there's nothing special we have to do here.
	return 0;
}
