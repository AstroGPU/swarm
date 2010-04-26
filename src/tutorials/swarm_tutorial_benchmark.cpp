#include "swarm.h"
#include "swarmlog.h"
#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <boost/program_options.hpp>

#define PRINT_OUTPUT 0
#define PARANOID_ENERGY_CHECK 1


// Declare functions to demonstrate setting/accessing system state
void set_initial_conditions_for_demo(swarm::ensemble& ens);
#if PRINT_OUTPUT
void print_selected_systems_for_demo(swarm::ensemble& ens);
#endif

int main(int argc,  char **argv)
{
  using namespace swarm;
  namespace po = boost::program_options;

  // performance stopwatches
  stopwatch swatch_kernel_gpu, swatch_upload_gpu, swatch_download_gpu, swatch_temps_gpu, swatch_init_gpu;
  stopwatch swatch_kernel_cpu, swatch_upload_cpu, swatch_temps_cpu, swatch_init_cpu;
  stopwatch swatch_all;

  swatch_all.start(); // Start timer for entire program
  srand(42u);    // Seed random number generator, so output is reproducible

  // Parse command line arguements (making it easy to compare)
  po::options_description desc(std::string("Usage: ") + argv[0] + " \nOptions");
  desc.add_options()
    ("help,h", "produce help message")
    ("systems,s", po::value<int>(), "number of systems [1,30720]")
    ("num_bodies,n", po::value<int>(),  "number of bodies per system [3,10]")
    ("time,t", po::value<double>(), "time to integrate (0,62832.)")
    ("blocksize,b", po::value<int>(),  "number of threads per block {16,32,48,64,128}")
    ("precision,p", po::value<int>(),  "precision (1=double, 2=single, 3=mixed)")
    ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  //
  bool valid = true;
  std::cerr << "# Set integrator parameters (hardcoded in this demo, except for command line arguments.).\n";
  config cfg;
  cfg["integrator"] = "gpu_hermite"; // integrator name
  cfg["runon"] = "gpu";              // whether to run on cpu or gpu (must match integrator)
  cfg["time step"] = "0.0005";       // time step

  //  cfg["precision"] = "1";            // use double precision

  // Get values for config hashmap from command line arguements (or use defaults)
  {
  std::ostringstream precision_stream;
  if(vm.count("precision")) 
    {
      int prec = vm["precision"].as<int>();
      precision_stream <<  prec;
      if(!((prec==1)||(prec==2)||(prec==3)))
	 valid =false;
    }
  else
    precision_stream << 1; 
  cfg["precision"] = precision_stream.str();
  }
  {
    std::ostringstream blocksize_stream;
    if(vm.count("blocksize")) 
      {
	int bs = vm["blocksize"].as<int>();
	blocksize_stream << bs;
	if((bs<16)||(bs>128)||(bs%16!=0))
	 valid =false;
      }
    else
      blocksize_stream << 64; 
    cfg["threads per block"] = blocksize_stream.str();
  }
  
  // Get simple values from command line arguements (or use defaults)
  int nsystems = (vm.count("systems")) ? vm["systems"].as<int>() : 960;
  int nbodyspersystem = (vm.count("num_bodies")) ? vm["num_bodies"].as<int>() : 3;
  double dT = (vm.count("time")) ? vm["time"].as<double>() : 2.*M_PI;

  // Print help message if used inappropriately
  if (vm.count("help")||!(nsystems>=1)||!((nbodyspersystem>=3)&&(nbodyspersystem<=10))||!((dT>0.)&&(dT<=2.*M_PI*10000.+1.))||!valid) { std::cout << desc << "\n"; return 1; }

  // Now that we've retreived parameters, start timers for initialization
  swatch_init_cpu.start();
  swatch_init_gpu.start();
  std::cerr << "# Initialize ensemble on host to be used with GPU integration.\n";
  cpu_ensemble ens(nsystems, nbodyspersystem);
  
  std::cerr << "# Set integration duration for all systems.\n";
  ens.set_time_end_all(dT);
  ens.set_time_output_all(1, 1.01*dT);	// time of next output is after integration ends

  std::cerr << "# Parameters: systems= " << nsystems << " num_bodies= " << nbodyspersystem << " time= " << dT << " blocksize= " << cfg["threads per block"] << ".\n";
  std:: cerr << "# Initialize the library\n";
  swarm::init(cfg);
  swatch_init_cpu.stop();

  std:: cerr << "# Initialize the GPU integrator\n";
  std::auto_ptr<integrator> integ_gpu(integrator::create(cfg));
  cudaThreadSynchronize();  // Block until CUDA call completes
  swatch_init_gpu.stop();

  std::cerr << "# Set initial conditions.\n";
  set_initial_conditions_for_demo(ens);
  
#if PARANOID_ENERGY_CHECK
  // Calculate energy at beginning of integration
  std::vector<double> energy_init(ens.nsys());
  ens.calc_total_energy(&energy_init[0]);
#endif

#if PRINT_OUTPUT
  // Print initial conditions on CPU for use w/ GPU
  std::cerr << "# Print selected initial conditions for GPU.\n";
  print_selected_systems_for_demo(ens);
#endif
  
  std::cerr << "# Create identical ensemble on host for comparison w/ CPU.\n";
  swatch_upload_cpu.start();
  cpu_ensemble ens_check(ens);
  swatch_upload_cpu.stop();

#if PRINT_OUTPUT
  // Print initial conditions for checking w/ CPU 
  std::cerr << "# Print selected initial conditions for CPU.\n";
  print_selected_systems_for_demo(ens_check);	
#endif
  
  // Perform the integration on gpu
  std::cerr << "# Upload data to GPU.\n";
  cudaThreadSynchronize();  // Block until CUDA call completes
  swatch_upload_gpu.start();
  gpu_ensemble gpu_ens(ens);
  cudaThreadSynchronize();  // Block until CUDA call completes
  swatch_upload_gpu.stop();

  std::cerr << "# Integrate ensemble on GPU.\n";
  swatch_temps_gpu.start();
  integ_gpu->integrate(gpu_ens, 0.);				
  cudaThreadSynchronize();  // Block until CUDA call completes
  swatch_temps_gpu.stop();

  swatch_kernel_gpu.start();
  integ_gpu->integrate(gpu_ens, dT);				
  cudaThreadSynchronize();  // Block until CUDA call completes
  swatch_kernel_gpu.stop();
  std::cerr << "# GPU integration complete.\n";

  std::cerr << "# Download data to host.\n";
  swatch_download_gpu.start();
  ens.copy_from(gpu_ens);					
  cudaThreadSynchronize();  // Block until CUDA call completes
  swatch_download_gpu.stop();
  std::cerr << "# Download complete.\n";
  
  // Perform the integration on the cpu
  std:: cerr << "# Initialize the CPU integrator\n";
  cfg["integrator"] = "cpu_hermite";
  swatch_init_cpu.start();
  std::auto_ptr<integrator> integ_cpu(integrator::create(cfg));
  swatch_init_cpu.stop();

  std::cerr << "# Integrate a copy of ensemble on CPU for comparison.\n";
  swatch_temps_cpu.start();
  integ_cpu->integrate(ens_check, 0.);				
  swatch_temps_cpu.stop();
  swatch_kernel_cpu.start();
  integ_cpu->integrate(ens_check, dT);				
  swatch_kernel_cpu.stop();
  std::cerr << "# CPU integration complete.\n";
  
  swatch_all.stop();

#if PRINT_OUTPUT
  // Print results
  std::cerr << "# Print selected results from GPU's calculation.\n";
  print_selected_systems_for_demo(ens);
  std::cerr << "# Print selected results from CPU's calculation.\n";
  print_selected_systems_for_demo(ens_check);
#endif
  
#if PARANOID_ENERGY_CHECK
  // Check Energy conservation
  std::vector<double> energy_gpu_final(ens.nsys()), energy_cpu_final(ens.nsys());;
  ens.calc_total_energy(&energy_gpu_final[0]);
  ens_check.calc_total_energy(&energy_cpu_final[0]);
  double max_deltaE = 0., max_deltaE_gpu = 0., max_deltaE_cpu = 0.;
  for(int sysid=0;sysid<ens.nsys();++sysid)
    {
      double deltaE_gpu = (energy_gpu_final[sysid]-energy_init[sysid])/energy_init[sysid];
      double deltaE_cpu = (energy_cpu_final[sysid]-energy_init[sysid])/energy_init[sysid];
      double deltaE = std::max(fabs(deltaE_gpu),fabs(deltaE_cpu));
      if(deltaE_gpu>max_deltaE_gpu)
	{ max_deltaE_gpu = deltaE_gpu; }
      if(deltaE_cpu>max_deltaE_cpu)
	{ max_deltaE_cpu = deltaE_cpu; }
      if(deltaE>max_deltaE)
	{ max_deltaE = deltaE; }
      if(fabs(deltaE)>0.00001)
	std::cout << "# Warning: " << sysid << " dE/E (gpu)= " << deltaE_gpu << " dE/E (cpu)= " << deltaE_cpu << '\n';
'\n';
    }
  std::cout.flush();
  std::cerr << "# Max dE/E (gpu)= " << max_deltaE_gpu << "  Max dE/E (cpu)= " << max_deltaE_cpu << "\n";
#endif  


  std::cerr << "# Time (all, combined): " << swatch_all.getTime()*1000. << " ms.\n";
  std::cerr << "# Time (init lib)     : " << swatch_init_gpu.getTime()*1000. << " ms on GPU,   " << swatch_init_cpu.getTime()*1000. << " ms on CPU.\n";
  std::cerr << "# Time (upload)       : " << swatch_upload_gpu.getTime()*1000. << " ms on GPU,   " << swatch_upload_cpu.getTime()*1000. << " ms on CPU.\n";
  std::cerr << "# Time (0th step)     : " << swatch_temps_gpu.getTime()*1000. << " ms on GPU,   " << swatch_temps_cpu.getTime()*1000. << " ms on CPU.\n";
  std::cerr << "# Time (integration)  : " << swatch_kernel_gpu.getTime()*1000. << " ms on GPU,   " << swatch_kernel_cpu.getTime()*1000. << " ms on CPU.\n";
  std::cerr << "# Time (download)     : " << swatch_download_gpu.getTime()*1000. << " ms on GPU,   " << 0. << " ms on CPU.\n";

  
  std::cerr << "# Speed up (kernel only)         : " << swatch_kernel_cpu.getTime()/swatch_kernel_gpu.getTime() << ".\n";
  std::cerr << "# Speed up (w/ mem transfer)     : " << (swatch_upload_cpu.getTime()+swatch_temps_cpu.getTime()+swatch_kernel_cpu.getTime())/(swatch_upload_gpu.getTime()+swatch_temps_gpu.getTime()+swatch_kernel_gpu.getTime()+swatch_download_gpu.getTime()) << ".\n";

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
	  x  =  rmag*cos(theta); y  = rmag*sin(theta); z  = 0;
	  vx = -vmag*sin(theta); vy = vmag*cos(theta); vz = 0.;
	  
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
  std::cout.flush();
}

