#include "swarm.h"
#include <iostream>

void set_initial_conditions_for_demo(ensemble& ens);
void print_selected_systems_for_demo(ensemble& ens);


int main(int argc, const char **argv)
{
	std::cerr << "Set integrator parameters\n";
	config cfg;
	cfg["integrator"] = "gpu_hermite"; // integrator name
	cfg["h"] = "0.0005";                 // time step
	cfg["precision"] = "1";            // use double precision

	std:: cerr << "Initialize the GPU integrator\n";
	std::auto_ptr<integrator> integ_gpu(integrator::create(cfg));

	std::cerr << "Initialize ensemble on host to be used with GPU integration.\n";
	unsigned int nsystems = 100, nbodyspersystem = 3;
	cpu_ensemble ens(nsystems, nbodyspersystem);

	std::cerr << "Set initial conditions.\n";
	set_initial_conditions_for_demo(ens);

	// Print initial conditions on CPU for use w/ GPU
	std::cerr << "Print selected initial conditions for GPU.\n";
	print_selected_systems_for_demo(ens);

	std::cerr << "Create identical ensemble on host to check w/ CPU.\n";
	cpu_ensemble ens_check(ens);

	// Print initial conditions for checking w/ CPU 
	std::cerr << "Print selected initial conditions for CPU.\n";
	print_selected_systems_for_demo(ens_check);	

	std::cerr << "Set integration duration for all systems.\n";
	double dT = 1.*2.*M_PI;
	ens.set_time_end_all(dT);

	// Perform the integration on gpu
	std::cerr << "Upload data to GPU.\n";
	gpu_ensemble gpu_ens(ens);
	std::cerr << "Integrate ensemble on GPU.\n";
	integ_gpu->integrate(gpu_ens, dT);				
	std::cerr << "Download data to host.\n";
	ens.copy_from(gpu_ens);					
	std::cerr << "GPU integration complete.\n";

	// Perform the integration on the cpu
	std:: cerr << "Initialize the CPU integrator\n";
	cfg["integrator"] = "cpu_hermite";
	std::auto_ptr<integrator> integ_cpu(integrator::create(cfg));
	std::cerr << "Integrate a copy of ensemble on CPU to check.\n";
	integ_cpu->integrate(ens_check, dT);				
	std::cerr << "CPU integration complete.\n";

	// Print results
	std::cerr << "Print selected results from GPU's calculation.\n";
	print_selected_systems_for_demo(ens);
	std::cerr << "Print selected results from CPU's calculation.\n";
	print_selected_systems_for_demo(ens_check);
	
	// both the integrator & the ensembles are automatically deallocated on exit
	// so there's nothing special we have to do here.
	return 0;
}




void set_initial_conditions_for_demo(ensemble& ens) 
{
  std::cerr << "Set initial time for all systems.\n";
  ens.set_time_all(0.);	

  for(unsigned int sys=0;sys<ens.nsys();++sys)
    {
      // set sun to unit mass and at origin
      double mass_sun = 1.;
      double x=0, y=0, z=0, vx=0, vy=0, vz=0;
      ens.set_body(sys, 0, mass_sun, x, y, z, vx, vy, vz);
      
      // add near-Jupiter-mass planets on nearly circular orbits
      for(unsigned int bod=1;bod<ens.nbod();++bod)
	{
	  float mass_planet = 0.001;
	  double rmag = pow(1.4,bod-1);
	  double vmag = sqrt(mass_sun/rmag);
	  double theta = (2.*M_PI*rand())/static_cast<double>(RAND_MAX);
	  x  = rmag*cos(theta); y  = rmag*sin(theta); z  = 0;
	  vx = vmag*sin(theta); vy = vmag*cos(theta); vz = 0.;
	  
	  // assign body a mass, position and velocity
	  ens.set_body(sys, bod, mass_planet, x, y, z, vx, vy, vz);
	}
    }
}

void print_selected_systems_for_demo(ensemble& ens)
{
  std::streamsize cout_precision_old = std::cout.precision();
  std::cout.precision(10);
  unsigned int nprint = 1;
  for(unsigned int systemid = 0; systemid< nprint; ++systemid)
    {
      std::cout << "sys= " << systemid << " time= " << ens.time(systemid) << "\n";
      for(unsigned int bod=0;bod<ens.nbod();++bod)
	{
	  std::cout << "body= " << bod << ": pos= (" << ens.x(systemid, bod) << ", " <<  ens.y(systemid, bod) << ", " << ens.z(systemid, bod) << ") vel= (" << ens.vx(systemid, bod) << ", " <<  ens.vy(systemid, bod) << ", " << ens.vz(systemid, bod) << ").\n";
	}
    }
  std::cout.precision(cout_precision_old);
}

