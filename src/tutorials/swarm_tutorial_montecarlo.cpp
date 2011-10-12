/*************************************************************************
 * Copyright (C) 2009-2010 by Eric Ford & the Swarm-NG Development Team  *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 3 of the License.        *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program; if not, write to the                         *
 * Free Software Foundation, Inc.,                                       *
 * 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ************************************************************************/

/*! \file swarm_tutorial_montecarlo.cpp
 *  \brief program for Monte Carlo integrations based on user specified initial conditions
*/

#define ASTROCENTRIC 1
#define BARRYCENTRIC 0 // Not implemented
#define JACOBI 0       // Not implemented

#include "swarm/swarm.h"
//#include "swarmlog.h"

#include <iostream>
#include <sstream>
#include <vector>
#include <memory>

#include "random.h"
#include "swarm/kepler.h"

#define PARANOID_CPU_CHECK 0  // WARNING:  Setting to 1 takes a very long time

void set_initial_conditions_for_demo(swarm::ensemble& ens, const swarm::config& cfg);
void print_selected_systems_for_demo(swarm::ensemble& ens);


int main(int argc, const char **argv)
{
  using namespace swarm;
  srand(42u);    // Seed random number generator, so output is reproducible

  // Get name of configuration file
  if(argc != 2)
    {
      std::cerr << "Usage: " << argv[0] << " <integrator.cfg>\n";
      return -1;
    }
  std::string icfgfn = argv[1];

  config cfg;
  //  EBF: Replaced w/ Saleh's new syntax
  //  load_config(cfg,icfgfn);
  cfg = config::load(icfgfn, cfg);
	       
   // Initialize Swarm library. Basically, it initializes CUDA system and default logging system.
  std:: cerr << "Initialize the library\n";
  swarm::init(cfg);

  std:: cerr << "Initialize the GPU integrator\n";
  // EBF: Replaced w/ Saleh's new syntax
  //  std::auto_ptr<integrator> integ_gpu(integrator::create(cfg));
  // Select and create the integrator. While you can create an integrator by calling its constructor directly.
  // It is recommended to use integrator::create since it gives you more flexibility at runtime.
  // The create function looks in the swarm library and finds the integrator you requested from 
  // the list of available plugins. 
  Pintegrator integ = integrator::create(cfg);
  
  std::cerr << "# Create ensemble on host to be used with CPU integration.\n";
  unsigned int nsystems = cfg.optional("num_systems",1024);
  unsigned int nbodyspersystem = cfg.optional("num_bodies",3);
  // EBF: Replace'd w/ Saleh's new syntax
  // cpu_ensemble ens(nsystems, nbodyspersystem);
  // defaultEnsemble ens = generate_ensemble(cfg);
  hostEnsemble ens_check = hostEnsemble::create( nbodyspersystem, nsystems );
  // std::cerr << "Set initial conditions on CPU.\n";
  set_initial_conditions_for_demo(ens_check,cfg);
 
   std::vector<double> energy_init(nsystems), energy_final(nsystems);
  ens_check.calc_total_energy(&energy_init[0]);
  
#if 1 // PARANOID_CPU_CHECK
  std::cerr << "Create identical ensemble on host to check w/ CPU.\n";
  // EBF: Replaced w/ Saleh's new syntax
  //  cpu_ensemble ens_check(ens);
  defaultEnsemble ens = ens_check.clone();

  
  // Print initial conditions for checking w/ CPU 
  std::cerr << "Print selected initial conditions for CPU.\n";
  print_selected_systems_for_demo(ens_check);	
#endif
  
  // Print initial conditions on CPU for use w/ GPU
  std::cerr << "Print selected initial conditions for upload to GPU.\n";
  print_selected_systems_for_demo(ens);
  

  // Now we set-up the integrator for integrating.
  //
  // First set the ensemble. For a GPU integrator, the GPU memory will be allocated.
  integ->set_ensemble(ens); 

  std::cerr << "Set integration duration for all systems.\n";
  // EBF: This was moved to set_initial_conditions_for_demo
  // double Tinit = cfg.optional("time_init",0.);
  double destination_time = cfg.optional("integration end",10.);
  // EBF: Replaced w/ Saleh's new syntax
  //  ens.set_time_end_all(Tend);
  // Now set the destination time where we want to stop the integration.
  integ->set_destination_time ( destination_time );
  // Need to synchronize because \ref integrator::set_ensemble may upload data to GPU.
  SYNC;

  // EBF: This is moved to the monitor for log_time_interval
  //  double Toutputstep = cfg.optional("output interval",destination_time*1.01);
  //  ens.set_time_output_all(0,Tinit);  // Time of first output (immediate)
  //  ens.set_time_output_all(1,Toutputstep);  // output interval
  // Perform logging if needed
  //  swarm::log::output_systems_needing_output(hlog, ens);

  // EBF: Disabled as this is now done automatically
  // Perform the integration on gpu
  //  std::cerr << "Upload data to GPU.\n";
  //  gpu_ensemble gpu_ens(ens);
  
  std::cerr << "Integrate ensemble on GPU.\n";
  // EBF: Replaced w/ Saleh's new syntax
  /*
  //  HACK:  Loop so that buffer does not overflow for hermite
  int MaxEntriesBeforeHermiteBufferOverFlows = 4; // static_cast<int>(std::floor(150000/(ens.nsys()*ens.nbod())));
  int nSubSections = static_cast<int>(std::ceil((Tend-Tinit)/Toutputstep))/MaxEntriesBeforeHermiteBufferOverFlows+1;
  for(int i=1;i<=nSubSections;++i)
    {
    double tmpTend = Tinit+(Tend-Tinit)*i/nSubSections;
    integ_gpu->integrate(gpu_ens, tmpTend);				
    }
  */
  // Now that everything is set-up, it is safe to pull the trigger and 
  // call the integrate method on integrator. Note that since we didn't set
  // a logger for the integrator, it will use the default logging system.
  integ->integrate();
  // Need to synchronize because integrate is a GPU call.
  SYNC;

  // EBF: Disabled as this is now done automatically
  //  std::cerr << "Download data to host.\n";
  //  ens.copy_from(gpu_ens);					
  std::cerr << "GPU integration complete.\n";
  
#if PARANOID_CPU_CHECK
#error Not implemented yet
  // Perform the integration on the cpu
  std:: cerr << "Initialize the CPU integrator\n";
  cfg["integrator"] = "cpu_hermite";
  std::auto_ptr<integrator> integ_cpu(integrator::create(cfg));
  std::cerr << "Integrate a copy of ensemble on CPU to check.\n";
  integ_cpu->integrate(ens_check, Tend);				
  std::cerr << "CPU integration complete.\n";
#endif
  
  // Check Energy conservation
  ens.calc_total_energy(&energy_final[0]);
  double max_deltaE = 0;
  for(int sysid=0;sysid<ens.nsys();++sysid)
    {
      double deltaE = (energy_final[sysid]-energy_init[sysid])/energy_init[sysid];
      if(fabs(deltaE)>max_deltaE)
	{ max_deltaE = fabs(deltaE); }
      if(fabs(deltaE)>0.00001)
	std::cout << "# Warning: " << sysid << " dE/E= " << deltaE << '\n';
    }
  std::cerr << "# Max dE/E= " << max_deltaE << "\n";
  
  // Print results
  std::cerr << "Print selected results from GPU's calculation.\n";
  print_selected_systems_for_demo(ens);
#if PARANOID_CPU_CHECK
  std::cerr << "Print selected results from CPU's calculation.\n";
  print_selected_systems_for_demo(ens_check);
#endif
  // both the integrator & the ensembles are automatically deallocated on exit
  // so there's nothing special we have to do here.
  return 0;
}




void set_initial_conditions_for_demo(swarm::ensemble& ens, const swarm::config& cfg) 
{
  using namespace swarm;
  std::cerr << "Set initial time for all systems = ";
  double time_init = cfg.optional("time_init",0.);
  for(unsigned int sys=0;sys<ens.nsys();++sys)
    {
      ens.set_time(sys,time_init);	
    }
  std::cerr << time_init << ".\n";

  bool use_jacobi = cfg.optional("use_jacobi",0);

  for(unsigned int sys=0;sys<ens.nsys();++sys)
    {
      // set sun to unit mass and at origin
      float mass_star = cfg.optional("mass_star",1.);
      double x=0, y=0, z=0, vx=0, vy=0, vz=0;
      ens.set_body(sys, 0, mass_star, x, y, z, vx, vy, vz);
      
	  double x_t = ens.x(sys,0);
	  double y_t = ens.y(sys,0);
	  double z_t = ens.z(sys,0);
	  double vx_t = ens.vx(sys,0);
	  double vy_t = ens.vy(sys,0);
	  double vz_t = ens.vz(sys,0);

	  std::cout << " x= " << x << "=" << x_t << " ";
	  std::cout << " y= " << y << "=" << y_t << " ";
	  std::cout << " z= " << z << "=" << z_t << " ";
	  std::cout << "vx= " << vx << "=" << vx_t << " ";
	  std::cout << "vy= " << vy << "=" << vy_t << " ";
	  std::cout << "vz= " << vz << "=" << vz_t << "\n";
			       
      double mass_enclosed = mass_star;
      for(unsigned int bod=1;bod<ens.nbod();++bod)
	{
	  float mass_planet = draw_value_from_config(cfg,"mass",bod,0.,100.);
	  mass_enclosed += mass_planet;

	  double a = draw_value_from_config(cfg,"a",bod,0.001,10000.);
	  double e = draw_value_from_config(cfg,"ecc",bod,0.,1.);
	  double i = draw_value_from_config(cfg,"inc",bod,-180.,180.);
	  double O = draw_value_from_config(cfg,"node",bod,-720.,720.);
	  double w = draw_value_from_config(cfg,"omega",bod,-720.,720.);
	  double M = draw_value_from_config(cfg,"meananom",bod,-720.,720.);
	  if(sys<10)
	    std::cout << "# Drawing sys= " << sys << " bod= " << bod << " m= " << mass_planet << "  a= " << a << " e= " << e << " i= " << i << " Omega= " << O << " omega= " << w << " M= " << M << '\n';

	  i *= M_PI/180.;
	  O *= M_PI/180.;
	  w *= M_PI/180.;
	  M *= M_PI/180.;

	  //	  double mass = use_jacobi ? mass_enclosed : mass_star+mass_planet;
	  double mass = mass_star+mass_planet;
	  calc_cartesian_for_ellipse(x,y,z,vx,vy,vz, a, e, i, O, w, M, mass);

	  double bx=0., by=0., bz=0., bvx=0., bvy=0., bvz=0.;
#if JACOBI
	  ens.get_barycenter(sys,bx,by,bz,bvx,bvy,bvz,bod-1);
#endif
	  x  += bx;	  y  += by;	  z  += bz;
	  vx += bvx;	  vy += bvy;	  vz += bvz;
	  
	  // assign body a mass, position and velocity
	  ens.set_body(sys, bod, mass_planet, x, y, z, vx, vy, vz);
	  x_t = ens.x(sys,bod);
	  y_t = ens.y(sys,bod);
	  z_t = ens.z(sys,bod);
	  vx_t = ens.vx(sys,bod);
	  vy_t = ens.vy(sys,bod);
	  vz_t = ens.vz(sys,bod);

	  std::cout << " x= " << x << "=" << x_t << " ";
	  std::cout << " y= " << y << "=" << y_t << " ";
	  std::cout << " z= " << z << "=" << z_t << " ";
	  std::cout << "vx= " << vx << "=" << vx_t << " ";
	  std::cout << "vy= " << vy << "=" << vy_t << " ";
	  std::cout << "vz= " << vz << "=" << vz_t << "\n";
	}  // end loop over bodies

      // Shift into barycentric frame
      ens.get_barycenter(sys,x,y,z,vx,vy,vz);
      for(unsigned int bod=0;bod<ens.nbod();++bod)
	{
	  ens.set_body(sys, bod, ens.mass(sys,bod), 
		       ens.x(sys,bod)-x, ens.y(sys,bod)-y, ens.z(sys,bod)-z, 
		       ens.vx(sys,bod)-vx, ens.vy(sys,bod)-vy, ens.vz(sys,bod)-vz);	  
	}  // end loop over bodies
    } // end loop over systems
}

void print_selected_systems_for_demo(swarm::ensemble& ens)
{
  using namespace swarm;
  std::streamsize cout_precision_old = std::cout.precision();
  std::cout.precision(10);
  unsigned int nprint = std::min(10,ens.nsys());
  for(unsigned int systemid = 0; systemid< nprint; ++systemid)
    {
      std::cout << "sys= " << systemid << " time= " << ens.time(systemid) << "\n";
      // EBF: Changed from float to double for Saleh's new data structure
      double star_mass = ens.mass(systemid,0);

      double mass_effective = star_mass;
      double bx, by, bz, bvx, bvy, bvz;
#if JACOBI
      ens.get_barycenter(systemid,bx,by,bz,bvx,bvy,bvz,bod-1);
#else
#if BARRYCENTRIC
      ens.get_barycenter(systemid,bx,by,bz,bvx,bvy,bvz);
#else
#if ASTROCENTRIC
      ens.get_body(systemid,0,star_mass,bx,by,bz,bvx,bvy,bvz);
#else
#error  "Must specify some reference frame!"
#endif
#endif
#endif

      {
	unsigned int bod = 0;
	  float mass = ens.mass(systemid,bod);
	  double x = ens.x(systemid,bod);
	  double y = ens.y(systemid,bod);
	  double z = ens.z(systemid,bod);
	  double vx = ens.vx(systemid,bod);
	  double vy = ens.vy(systemid,bod);
	  double vz = ens.vz(systemid,bod);

	  std::cout << "body= " << bod << ": mass= " << mass << " pos= " << x << ' ' << y << ' ' << z << " vel= " << vx << ' ' << vy << ' ' << vz << "\n";
     }

      for(unsigned int bod=1;bod<ens.nbod();++bod) // Skip star since printing orbits
	{
	  //	  std::cout << "body= " << bod << ": ";
	  //	  std::cout << "pos= (" << ens.x(systemid, bod) << ", " <<  ens.y(systemid, bod) << ", " << ens.z(systemid, bod) << ") vel= (" << ens.vx(systemid, bod) << ", " <<  ens.vy(systemid, bod) << ", " << ens.vz(systemid, bod) << ").\n";
	  float mass = ens.mass(systemid,bod);
#if JACOBI
	  ens.get_barycenter(systemid,bx,by,bz,bvx,bvy,bvz,bod-1);
	  mass_effective += mass;
#else
#if BARRYCENTRIC
	  //	  ens.get_barycenter(systemid,bx,by,bz,bvx,bvy,bvz);
	  mass_effective = star_mass + mass;
#else
#if ASTROCENTRIC
	  //	  ens.get_body(systemid,0,star_mass,bx,by,bz,bvx,bvy,bvz);
	  mass_effective = star_mass + mass;
#else
#error    "Must specify some reference frame!"
#endif
#endif
#endif
	      double x = ens.x(systemid,bod)-bx;
	      double y = ens.y(systemid,bod)-by;
	      double z = ens.z(systemid,bod)-bz;
	      double vx = ens.vx(systemid,bod)-bvx;
	      double vy = ens.vy(systemid,bod)-bvy;
	      double vz = ens.vz(systemid,bod)-bvz;

	      std::cout << " body= " << bod << ": m= " << mass << " pos= " << x << ' ' << y << ' ' << z << " vel= " << vx << ' ' << vy << ' ' << vz << "\n";
	      double a, e, i, O, w, M;
	      calc_keplerian_for_cartesian(a,e,i,O,w,M, x,y,z,vx,vy,vz, mass_effective);
	      i *= 180/M_PI;
	      O *= 180/M_PI;
	      w *= 180/M_PI;
	      M *= 180/M_PI;
	      std::cout << "body= " << bod << ": m= " << mass << " a= " << a << " e= " << e << " i= " << i << " Omega= " << O << " omega= " << w << " M= " << M << "\n";


	}
    }
  std::cout.precision(cout_precision_old);
  std::cout << std::flush;
}

