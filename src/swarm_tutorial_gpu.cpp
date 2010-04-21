#include "swarm.h"
#include "swarmlog.h"
#include <iostream>
#include <memory>

// Declare functions to demonstrate setting/accessing system state
void set_initial_conditions_for_demo(swarm::ensemble& ens);
void print_selected_systems_for_demo(swarm::ensemble& ens);

int main(int argc, const char **argv)
{
  using namespace swarm;
  std::cerr << "Set integrator parameters (hardcoded in this demo).\n";
  config cfg;
  cfg["integrator"] = "gpu_rk4"; // integrator name
  cfg["h"] = "0.0005";               // time step
  cfg["precision"] = "1";            // use double precision
  cfg["rmax"] = "1000.";              // count the planet as "ejected" if it ventures beyond this radius (not all integrators support this)
  cfg["output"] = "null";            // store no output

  std:: cerr << "Initialize the library\n";
  swarm::init(cfg);

  std:: cerr << "Initialize the integrator\n";
  std::auto_ptr<integrator> integ(integrator::create(cfg));
  
  std::cerr << "Initialize ensemble on host to be used with GPU integration.\n";
  unsigned int nsystems = 128, nbodyspersystem = 3;
  cpu_ensemble ens(nsystems, nbodyspersystem);
  
  std::cerr << "Set initial conditions.\n";
  set_initial_conditions_for_demo(ens);
  
  std::cerr << "Print selected initial conditions for GPU.\n";
  print_selected_systems_for_demo(ens);
  
  std::cerr << "Set integration duration for all systems.\n";
  double dT = 1.*2.*M_PI;
  ens.set_time_end_all(dT);
  ens.set_time_output_all(1, 1.01*dT);	// time of next output is after integration ends -- effectively disable the outputs (not all integrators support this)
  
  std::cerr << "Upload data to GPU.\n";
  gpu_ensemble gpu_ens(ens);
  std::cerr << "Upload complete.\n";

  std::cerr << "Integrate ensemble on GPU.\n";
  integ->integrate(gpu_ens, dT);				
  std::cerr << "Integration complete.\n";

  std::cerr << "Download data to host.\n";
  ens.copy_from(gpu_ens);					
  std::cerr << "Download complete.\n";
  
  std::cerr << "Print selected results from GPU's calculation.\n";
  print_selected_systems_for_demo(ens);
  
  // both the integrator & the ensembles are automatically deallocated on exit
  // so there's nothing special we have to do here.
  return 0;
}



// Demonstrates how to assign initial conditions to a swarm::ensemble object
void set_initial_conditions_for_demo(swarm::ensemble& ens) 
{
  using namespace swarm;

  ens.set_time_all(0.);	  // Set initial time for all systems.

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
	  x  = rmag*cos(theta); y  = rmag*sin(theta); z  = 0;
	  vx = vmag*sin(theta); vy = vmag*cos(theta); vz = 0.;
	  
	  // assign body a mass, position and velocity
	  ens.set_body(sys, bod, mass_planet, x, y, z, vx, vy, vz);
	}
    }
}

// Demonstrates how to extract the state of each body from a swarm::ensemble object
void print_selected_systems_for_demo(swarm::ensemble& ens)
{
  using namespace swarm;
  std::streamsize cout_precision_old = std::cout.precision();    // Restore prcission of cout before calling this function
  std::cout.precision(10);  // Print at higher precission
  unsigned int nprint = 1;  // Limit output to first nprint system(s)
  for(unsigned int systemid = 0; systemid< nprint; ++systemid)
    {
      std::cout << "sys= " << systemid << " time= " << ens.time(systemid) << " nsteps= " << ens.nstep(systemid) << "\n";
      for(unsigned int bod=0;bod<ens.nbod();++bod)
	{
	  std::cout << "body= " << bod << ": mass= " << ens.mass(systemid, bod) << " position= (" << ens.x(systemid, bod) << ", " <<  ens.y(systemid, bod) << ", " << ens.z(systemid, bod) << ") velocity= (" << ens.vx(systemid, bod) << ", " <<  ens.vy(systemid, bod) << ", " << ens.vz(systemid, bod) << ").\n";
	}
    }
  std::cout.precision(cout_precision_old);  // Restore old precission to cout
}

