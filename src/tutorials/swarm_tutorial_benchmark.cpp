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

  // Set integrator parameters (no config file for this tutorial)
  // integrator, runon and time step are hard coded, rest can come from command line
  config cfg;
  cfg["integrator"] = "gpu_hermite"; // Set to use a GPU integrator
  cfg["runon"]      = "gpu";         // Set to runon GPU
  cfg["time step"] = "0.0005";       // time step
  bool valid = true;                 // Indicates whether cfg parameters are valid

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

  // Check that parameters from command line are ok
  if(!(nsystems>=1)||!(nsystems<=32720)) valid = false;
  if(!(nbodyspersystem>=3)||!(nbodyspersystem<=10)) valid = false;
  if(!(dT>0.)||!(dT<=2.*M_PI*10000.+1.)) valid = false;

  // Print help message if requested or invalid parameters
  if (vm.count("help")||!valid) { std::cout << desc << "\n"; return 1; }

  // Print parameters for this set of benchmarks
  std::cerr << "# Parameters: systems= " << nsystems << " num_bodies= " << nbodyspersystem << " time= " << dT << " precision= " << cfg["precision"] << " blocksize= " << cfg["threads per block"] << ".\n";

  // Start CPU & GPU timers for initialization
  swatch_init_cpu.start();
  swatch_init_gpu.start();

  // Initialize ensemble on host to be used with GPU integration.
  cpu_ensemble ens(nsystems, nbodyspersystem);
  
  ens.set_time_end_all(dT);  // Set integration duration for all systems.
  ens.set_time_output_all(1, 1.01*dT);	// Set time of next output to be after integration ends

  swarm::init(cfg);         // Initialize the library
  swatch_init_cpu.stop();   // Stop timer for cpu initialization

  // Initialize the GPU integrator
  std::auto_ptr<integrator> integ_gpu(integrator::create(cfg));
  cudaThreadSynchronize();  // Block until CUDA call completes
  swatch_init_gpu.stop();   // Stop timer for cpu initialization

  std::cerr << "# Set initial conditions.\n";
  set_initial_conditions_for_demo(ens);
  
#if PARANOID_ENERGY_CHECK
  // Calculate energy at beginning of integration
  std::vector<double> energy_init(ens.nsys());
  ens.calc_total_energy(&energy_init[0]);
#endif

#if PRINT_OUTPUT
  std::cerr << "# Print selected initial conditions for GPU.\n";
  print_selected_systems_for_demo(ens);
#endif
  
  swatch_upload_cpu.start();   // Start timer for copying initial conditions into new CPU ensemble
  cpu_ensemble ens_check(ens); // Make a copy of the CPU ensemble for comparison
  swatch_upload_cpu.stop();    // Stop timer for copying ICs into new CPU ensemble

#if PRINT_OUTPUT
  std::cerr << "# Print selected initial conditions for CPU.\n";
  print_selected_systems_for_demo(ens_check);	
#endif
  
  std::cerr << "# Upload data to GPU.\n";
  cudaThreadSynchronize();   // Block until CUDA call completes
  swatch_upload_gpu.start(); // Start timer for copyg initial conditions to GPU
  gpu_ensemble gpu_ens(ens); // Initialize GPU ensemble, incl. copying data from CPU
  cudaThreadSynchronize();   // Block until CUDA call completes
  swatch_upload_gpu.stop();  // Stop timer for copyg initial conditions to GPU

  std::cerr << "# Integrate ensemble on GPU.\n";
  swatch_temps_gpu.start();  // Start timer for 0th step on GPU
  integ_gpu->integrate(gpu_ens, 0.);  // a 0th step of dT=0 results in initialization of the integrator only
  cudaThreadSynchronize();   // Block until CUDA call completes
  swatch_temps_gpu.stop();   // Stop timer for 0th step on GPU

  swatch_kernel_gpu.start(); // Start timer for GPU integration kernel
  integ_gpu->integrate(gpu_ens, dT);  // Actually do the integration w/ GPU!			
  cudaThreadSynchronize();  // Block until CUDA call completes
  swatch_kernel_gpu.stop(); // Stop timer for GPU integration kernel  
  std::cerr << "# GPU integration complete.\n";

  std::cerr << "# Download data to host.\n";
  swatch_download_gpu.start();  // Start timer for downloading data from GPU
  ens.copy_from(gpu_ens);	// Download data from GPU to CPU		
  cudaThreadSynchronize();      // Block until CUDA call completes
  swatch_download_gpu.stop();   // Stop timer for downloading data from GPU
  std::cerr << "# Download complete.\n";
  
  // Get ready to perform the same integration on the cpu
  cfg["integrator"] = "cpu_hermite"; // change to CPU version of integrator
  cfg["runon"]      = "cpu";         // change to runon CPU
  swatch_init_cpu.start();           // restart timer for initializing CPU integrator
  // Initialize the CPU integrator
  std::auto_ptr<integrator> integ_cpu(integrator::create(cfg));
  swatch_init_cpu.stop();            // Stop timer for initializing CPU integrator

  std::cerr << "# Integrate a copy of ensemble on CPU for comparison.\n";
  swatch_temps_cpu.start();      // Start timer for 0th step on CPU  
  integ_cpu->integrate(ens_check, 0.);  // a 0th step of dT=0 results in initialization of the integrator only				
  swatch_temps_cpu.stop();       // Stop timer for 0th step on CPU

  swatch_kernel_cpu.start();     // Start timer for CPU integration kernel
  integ_cpu->integrate(ens_check, dT);	 // Actually do the integration w/ CPU!
  swatch_kernel_cpu.stop();      // Stop timer for CPU integration kernel

  swatch_all.stop();             // Stop timer for all calculations
  std::cerr << "# CPU integration complete.\n";

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

  std::cerr << "# Benchmark results:\n";
  std::cerr << "# Time (all, combined): " << swatch_all.getTime()*1000. << " ms.\n";
  std::cerr << "# Time (init lib)     : " << swatch_init_gpu.getTime()*1000. << " ms on GPU,   " << swatch_init_cpu.getTime()*1000. << " ms on CPU.\n";
  std::cerr << "# Time (upload)       : " << swatch_upload_gpu.getTime()*1000. << " ms on GPU,   " << swatch_upload_cpu.getTime()*1000. << " ms on CPU.\n";
  std::cerr << "# Time (0th step)     : " << swatch_temps_gpu.getTime()*1000. << " ms on GPU,   " << swatch_temps_cpu.getTime()*1000. << " ms on CPU.\n";
  std::cerr << "# Time (integration)  : " << swatch_kernel_gpu.getTime()*1000. << " ms on GPU,   " << swatch_kernel_cpu.getTime()*1000. << " ms on CPU.\n";
  std::cerr << "# Time (download)     : " << swatch_download_gpu.getTime()*1000. << " ms on GPU,   " << 0. << " ms on CPU.\n";

  
  std::cerr << "# Speed up (kernel only)         : " << swatch_kernel_cpu.getTime()/swatch_kernel_gpu.getTime() << ".\n";
  std::cerr << "# Speed up (w/ mem transfer)     : " << (swatch_upload_cpu.getTime()+swatch_temps_cpu.getTime()+swatch_kernel_cpu.getTime())/(swatch_upload_gpu.getTime()+swatch_temps_gpu.getTime()+swatch_kernel_gpu.getTime()+swatch_download_gpu.getTime()) << ".\n";

  std::cout << swatch_kernel_cpu.getTime()/swatch_kernel_gpu.getTime() << "\n";

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

