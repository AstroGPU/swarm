/*************
 *  Author : Saleh Dindar
 *
 *
 */
#include "utils.hpp"
#include <algorithm>
using std::max;
using namespace swarm;

void generate_ensemble(config& cfg, cpu_ensemble& ens)  {

	double duration = atof(cfg["duration"].c_str());
	int nsys = atoi(cfg["nsys"].c_str());
	int nbod = atoi(cfg["nbod"].c_str());


	ens.reset(nsys,nbod,true);

	ens.set_time_all(0.);	  // Set initial time for all systems.
	ens.set_time_end_all(duration);  // Set integration duration for all systems.
	ens.set_time_output_all(1, 1.01*duration);	// Set time of next output to be after integration ends

	for(unsigned int sys=0;sys<ens.nsys();++sys)
	{
		// set sun to unit mass and at origin
		double mass_sun = 1.;
		double x=0, y=0, z=0, vx=0, vy=0, vz=0;
		ens.set_body(sys, 0, mass_sun, x, y, z, vx, vy, vz);

		// add near-Jupiter-mass planets on nearly circular orbits
		for(unsigned int bod=1;bod<ens.nbod();++bod)
		{
			float mass_planet = 0.001; // approximately (mass of Jupiter)/(mass of sun)
			double rmag = pow(1.4,bod-1);  // semi-major axes exceeding this spacing results in systems are stable for nbody=3 and mass_planet=0.001
			double vmag = sqrt(mass_sun/rmag);  // spped for uniform circular motion
			double theta = (2.*M_PI*rand())/static_cast<double>(RAND_MAX);  // randomize initial positions along ecah orbit
			x  =  rmag*cos(theta); y  = rmag*sin(theta); z  = 0;
			vx = -vmag*sin(theta); vy = vmag*cos(theta); vz = 0.;

			// assign body a mass, position and velocity
			ens.set_body(sys, bod, mass_planet, x, y, z, vx, vy, vz);
		}
	}
}

double find_max_energy_conservation_error(cpu_ensemble& ens, cpu_ensemble& reference_ensemble ) {
	std::vector<double> energy_init(reference_ensemble.nsys());
#if 1
	// Shift into center-of-mass frame
	for(unsigned int i=0; i<reference_ensemble.nsys() ; ++i)
	  {
	    double cx=0., cy=0., cz=0., cvx=0., cvy=0., cvz=0., msum=0.;
	    for(unsigned int j=0; j<reference_ensemble.nbod(); ++j)
	      {
		msum += reference_ensemble.mass(i,j);
		cx   += reference_ensemble.mass(i,j) * reference_ensemble.x(i,j);
		cy   += reference_ensemble.mass(i,j) * reference_ensemble.y(i,j);
		cz   += reference_ensemble.mass(i,j) * reference_ensemble.z(i,j);
		cvx   += reference_ensemble.mass(i,j) * reference_ensemble.vx(i,j);
		cvy   += reference_ensemble.mass(i,j) * reference_ensemble.vy(i,j);
		cvz   += reference_ensemble.mass(i,j) * reference_ensemble.vz(i,j);
	      }
	    cx /= msum; cy /= msum; cz /= msum; cvx /= msum; cvy /= msum; cvz /= msum;
	    for(unsigned int j=0; j<reference_ensemble.nbod(); ++j)
	      {
		reference_ensemble.x(i,j) -= cx;
		reference_ensemble.y(i,j) -= cy;
		reference_ensemble.z(i,j) -= cz;
		reference_ensemble.vx(i,j) -= cvx;
		reference_ensemble.vy(i,j) -= cvy;
		reference_ensemble.vz(i,j) -= cvz;
	      }
	  }
#endif
	reference_ensemble.calc_total_energy(&energy_init[0]);
	std::vector<double> energy_final(ens.nsys());
#if 1
	// Shift into center-of-mass frame
	for(unsigned int i=0; i<ens.nsys() ; ++i)
	  {
	    double cx=0., cy=0., cz=0., cvx=0., cvy=0., cvz=0., msum=0.;
	    for(unsigned int j=0; j<ens.nbod(); ++j)
	      {
		msum += ens.mass(i,j);
		cx   += ens.mass(i,j) * ens.x(i,j);
		cy   += ens.mass(i,j) * ens.y(i,j);
		cz   += ens.mass(i,j) * ens.z(i,j);
		cvx   += ens.mass(i,j) * ens.vx(i,j);
		cvy   += ens.mass(i,j) * ens.vy(i,j);
		cvz   += ens.mass(i,j) * ens.vz(i,j);
	      }
	    cx /= msum; cy /= msum; cz /= msum; cvx /= msum; cvy /= msum; cvz /= msum;
	    for(unsigned int j=0; j<ens.nbod(); ++j)
	      {
		ens.x(i,j) -= cx;
		ens.y(i,j) -= cy;
		ens.z(i,j) -= cz;
		ens.vx(i,j) -= cvx;
		ens.vy(i,j) -= cvy;
		ens.vz(i,j) -= cvz;
	      }
	  }
#endif
	ens.calc_total_energy(&energy_final[0]);
	double max_deltaE = 0.;
	for(int sysid=0;sysid<ens.nsys();++sysid)
	{
		double deltaE = (energy_final[sysid]-energy_init[sysid])/energy_init[sysid];
		max_deltaE = max(deltaE, max_deltaE);
	}
	return max_deltaE;
}

bool validate_configuration(config& cfg){
  bool valid = true;                 // Indicates whether cfg parameters are valid
  int nsystems = atoi(cfg["nsys"].c_str());
  int nbodyspersystem = atoi(cfg["nbod"].c_str());
  double dT = atof(cfg["duration"].c_str());
  int bs = atoi(cfg["blocksize"].c_str()) ;
  int prec = atoi(cfg["precision"].c_str()) ;


  // Check that parameters from command line are ok
  if(!((prec==1)||(prec==2)||(prec==3))) valid =false;
  if((bs<8)||(bs>512)) valid =false;
  if(!(nsystems>=1)||!(nsystems<=32720)) valid = false;
  if(!(nbodyspersystem>=3)||!(nbodyspersystem<=10)) valid = false;
  if(!(dT>0.)||!(dT<=2.*M_PI*1000000.+1.)) valid = false;

  return valid;
}
