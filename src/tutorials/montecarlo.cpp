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

/*! \file montecarlo.cpp
 *  \brief Implement Monte Carlo simulation to find planetary systems and generate ensemble. 
 *
 */


/** 
 * In writing this monte carlo simulation which is supposed to find planetary
 * systems using Monte Carlo simulations, I used the old monte carlo code.
 * It is useful, the tricky part is generating the ensemble.
 *
 *
 */
#include <iostream>
#include <fstream>
#include <math.h>
#include <signal.h>

#include "swarm/swarm.h"
#include "swarm/snapshot.hpp"
#include "swarm/log/log.hpp"
#include "random.hpp"
#include "kepler.hpp"

#define SYNC cudaThreadSynchronize()


using namespace swarm;
using namespace std;

config cfg;

void inspect(defaultEnsemble &ens, const int sys, const int bod ) 
{
	  fprintf(stderr,"%d %d: %lg (%lg %lg %lg) (%lg %lg %lg) \n", sys, bod,
			  ens[sys][bod].mass(),
			  ens[sys][bod][0].pos(),
			  ens[sys][bod][1].pos(),
			  ens[sys][bod][2].pos(),
			  ens[sys][bod][0].vel(),
			  ens[sys][bod][1].vel(),
			  ens[sys][bod][2].vel()
			  );
}

void generate_initial_conditions_for_system(const config& cfg, defaultEnsemble &ens, const int sysidx, const int sysid) 
{
  double time_init = cfg.optional("time_init", 0.0);
  const bool use_jacobi = cfg.optional("use_jacobi", 0);
  bool keplerian_ephemeris = cfg.optional("use_keplerian", 1);
  bool transit_ephemeris = cfg.optional("use_transit", 0);
  if(transit_ephemeris) keplerian_ephemeris = 0;
  assert(!(keplerian_ephemeris && transit_ephemeris));  // find a better way

  ens[sysidx].id() = sysid;
  ens[sysidx].time() = time_init;
  ens[sysidx].set_active();
  // set sun to unit mass and at origin
  double mass_star = cfg.optional("mass_star", 1.);
  double x=0, y=0, z=0, vx=0, vy=0, vz=0;
  ens.set_body(sysidx, 0, mass_star, x, y, z, vx, vy, vz);
  if(ens[sysidx][0].num_attributes()>=1)
    ens[sysidx][0].attribute(0) = cfg.optional("radius_star", 1.) * 0.00464912633;

  double mass_enclosed = mass_star;
  for(unsigned int bod=1;bod<ens.nbod();++bod)
    {
      double mass_planet = draw_value_from_config(cfg,"mass",bod,0.,mass_star);
      mass_enclosed += mass_planet;
      
      if(ens[sysidx][bod].num_attributes()>=1)
	ens[sysidx][bod].attribute(0) = draw_value_from_config(cfg,"radius",bod,0.,200.) * 4.26349283e-5;


      double a, e, i, O, w, M;
      if(transit_ephemeris)
	{
	  double period = draw_value_from_config(cfg,"period",bod,0.2,365250.);
	  double epoch = draw_value_from_config(cfg,"epoch",bod,0.2,365250.);
	  a = pow((period/365.25)*(period/365.25)*mass_enclosed,1.0/3.0);
	  e = draw_value_from_config(cfg,"ecc",bod,0.,1.);
	  i = draw_value_from_config(cfg,"inc",bod,-180.,180.);
	  O = draw_value_from_config(cfg,"node",bod,-720.,720.);
	  w = draw_value_from_config(cfg,"omega",bod,-720.,720.);
	  
	  i *= M_PI/180.;
	  O *= M_PI/180.;
	  w *= M_PI/180.;
	  
	  double T = (0.5*M_PI-w)+2.0*M_PI*((epoch/period)-floor(epoch/period));
	  double E = atan2(sqrt(1.-e)*(1.+e)*sin(T),e+cos(T));
	  M = E-e*sin(E);
	}
      else if (keplerian_ephemeris)
	{
	  a = draw_value_from_config(cfg,"a",bod,0.001,10000.);
	  e = draw_value_from_config(cfg,"ecc",bod,0.,1.);
	  i = draw_value_from_config(cfg,"inc",bod,-180.,180.);
	  O = draw_value_from_config(cfg,"node",bod,-720.,720.);
	  w = draw_value_from_config(cfg,"omega",bod,-720.,720.);
	  M = draw_value_from_config(cfg,"meananom",bod,-720.,720.);
	  
	  i *= M_PI/180.;
	  O *= M_PI/180.;
	  w *= M_PI/180.;
	  M *= M_PI/180.;
	}
      
      double mass = use_jacobi ? mass_enclosed : mass_star+mass_planet;
      if(cfg.count("verbose")&&(sysid<10))
	std::cout << "# Drawing sysid= " << sysid << " bod= " << bod << ' ' << mass_planet << "  " << a << ' ' << e << ' ' << i*180./M_PI << ' ' << O*180./M_PI << ' ' << w*180./M_PI << ' ' << M*180./M_PI << '\n';
      
      calc_cartesian_for_ellipse(x,y,z,vx,vy,vz, a, e, i, O, w, M, mass);
      
      //  printf("%d %d: %lg (%lg %lg %lg) (%lg %lg %lg) \n", sysidx, bod, mass_planet, x,y,z,vx,vy,vz);
      
      if(use_jacobi)
	{
	  double bx, by, bz, bvx, bvy, bvz;
	  ens.get_barycenter(sysidx,bx,by,bz,bvx,bvy,bvz,bod-1);
	  x  += bx;	  y  += by;	  z  += bz;
	  vx += bvx;	  vy += bvy;	  vz += bvz;
	}
      
      // assign body a mass, position and velocity
      ens.set_body(sysidx, bod, mass_planet, x, y, z, vx, vy, vz);
      
      if(cfg.count("verbose")&&(sysid<10))
	{
	  double x_t = ens.x(sysidx,bod);
	  double y_t = ens.y(sysidx,bod);
	  double z_t = ens.z(sysidx,bod);
	  double vx_t = ens.vx(sysidx,bod);
	  double vy_t = ens.vy(sysidx,bod);
	  double vz_t = ens.vz(sysidx,bod);
	  
	  std::cout << " x= " << x << "=" << x_t << " ";
	  std::cout << " y= " << y << "=" << y_t << " ";
	  std::cout << " z= " << z << "=" << z_t << " ";
	  std::cout << "vx= " << vx << "=" << vx_t << " ";
	  std::cout << "vy= " << vy << "=" << vy_t << " ";
	  std::cout << "vz= " << vz << "=" << vz_t << "\n";
	}
      
    }  // end loop over bodies
  
  // Shift into barycentric frame
  ens.get_barycenter(sysidx,x,y,z,vx,vy,vz);
  for(unsigned int bod=0;bod<ens.nbod();++bod)
    {
      ens.set_body(sysidx, bod, ens.mass(sysidx,bod), 
		   ens.x(sysidx,bod)-x, ens.y(sysidx,bod)-y, ens.z(sysidx,bod)-z, 
		   ens.vx(sysidx,bod)-vx, ens.vy(sysidx,bod)-vy, ens.vz(sysidx,bod)-vz);	  
    }  // end loop over bodies
}

/**
 * This procedure was copied from the other Monte Carlo simulation code
 *
 * I had to change all the floats to double because in the new swarm we
 * use double for masses. Also, our ensemble classes now use reference
 * counted pointers and we create them in a functional way than old 
 * fill-this-for-me-please pointers(references). The generate name
 * makes more sense that the set_initial_conditions. I only added
 * one line for generating the ensemble. In the other monte carlo, 
 * the ensemble is generated first and then passed to this function
 * to fill it in.
 *
 *
 */
defaultEnsemble generate_ensemble_with_randomized_initial_conditions(const config& cfg) 
{
  defaultEnsemble ens = defaultEnsemble::create( cfg.require("nbod",0), cfg.require("nsys",0) );
  std::cerr << "# nsystems = " << ens.nsys() << " nbodies/system = " << ens.nbod() << "\n";

  //  std::cerr << "# Set initial time for all systems = ";
  double time_init = cfg.optional("time_init", 0.0);
  //  std::cerr << time_init << ".\n";

  const bool use_jacobi = cfg.optional("use_jacobi", 0);
  bool keplerian_ephemeris = cfg.optional("use_keplerian", 1);
  bool transit_ephemeris = cfg.optional("use_transit", 0);
  if(transit_ephemeris) keplerian_ephemeris = 0;
  assert(!(keplerian_ephemeris && transit_ephemeris));  // find a better way

  for(unsigned int sys=0;sys<ens.nsys();++sys)
    {
      generate_initial_conditions_for_system(cfg,ens,sys,sys);
    } // end loop over systems
  return ens;
}


void print_system(const swarm::ensemble& ens, const int systemid, std::ostream &os = std::cout)
{
  enum {
    JACOBI, BARYCENTRIC, ASTROCENTRIC
  } COORDINATE_SYSTEM = BARYCENTRIC;
  const bool use_jacobi = cfg.optional("use_jacobi", 0);  
  if(use_jacobi) COORDINATE_SYSTEM = JACOBI;

  std::streamsize cout_precision_old = os.precision();
  os.precision(10);
  os << "sys_idx= " << systemid << " sys_id= " << ens[systemid].id() << " time= " << ens.time(systemid) << "\n";
      double star_mass = ens.mass(systemid,0);
      double mass_effective = star_mass;
      double bx, by, bz, bvx, bvy, bvz;
      switch(COORDINATE_SYSTEM) 
	{
	case JACOBI:
	  ens.get_barycenter(systemid,bx,by,bz,bvx,bvy,bvz,0);
	  break;
	  
	case BARYCENTRIC:
	  ens.get_barycenter(systemid,bx,by,bz,bvx,bvy,bvz);
	  break;
	  
	case ASTROCENTRIC:
	  ens.get_body(systemid,0,star_mass,bx,by,bz,bvx,bvy,bvz);
	  break;
	}
      
      for(unsigned int bod=1;bod<ens.nbod();++bod) // Skip star since printing orbits
	{
	  //	  std::cout << "pos= (" << ens.x(systemid, bod) << ", " <<  ens.y(systemid, bod) << ", " << ens.z(systemid, bod) << ") vel= (" << ens.vx(systemid, bod) << ", " <<  ens.vy(systemid, bod) << ", " << ens.vz(systemid, bod) << ").\n";
	  double mass = ens.mass(systemid,bod);
	  
	  switch(COORDINATE_SYSTEM) 
	    {
	    case JACOBI:
	      ens.get_barycenter(systemid,bx,by,bz,bvx,bvy,bvz,bod-1);
	      mass_effective += mass;
	      break;
	      
	    case BARYCENTRIC:
	      mass_effective = star_mass + mass;
	      break;
	      
	    case ASTROCENTRIC:
	      mass_effective = star_mass + mass;
	      break;
	    }
	  double x = ens.x(systemid,bod)-bx;
	  double y = ens.y(systemid,bod)-by;
	  double z = ens.z(systemid,bod)-bz;
	  double vx = ens.vx(systemid,bod)-bvx;
	  double vy = ens.vy(systemid,bod)-bvy;
	  double vz = ens.vz(systemid,bod)-bvz;
	  
	  double a, e, i, O, w, M;
	  calc_keplerian_for_cartesian(a,e,i,O,w,M, x,y,z,vx,vy,vz, mass_effective);
	  i *= 180/M_PI;
	  O *= 180/M_PI;
	  w *= 180/M_PI;
	  M *= 180/M_PI;
	  //	  os << " b= " << bod << " m= " << mass << " a= " << a << " e= " << e << " i= " << i << " Omega= " << O << " omega= " << w << " M= " << M << "\n";
	  os << ens[systemid].id() << " " << bod << " " << mass << " " << a << " " << e << " " << i << " " << O << " " << w << " " << M << "\n";
	}

  os.precision(cout_precision_old);
  os << std::flush;
}


void print_selected_systems(swarm::ensemble& ens, std::vector<unsigned int> systemindices, std::ostream &os = std::cout)
{
  for(unsigned int i=0; i<systemindices.size(); ++i)
      print_system(ens,systemindices[i], os);
}

void print_selected_systems_for_demo(swarm::ensemble& ens, unsigned int nprint, std::ostream &os = std::cout)
{
  if(nprint>ens.nsys()) nprint = ens.nsys();
  for(unsigned int systemid = 0; systemid< nprint; ++systemid)
      print_system(ens,systemid,os);
}


void write_stable_systems(defaultEnsemble &ens, defaultEnsemble &ens_init) 
{
  // find the stable ones and output the initial conditions for the stable
  // ones in keplerian coordinates
  //  std::cerr << "# Finding stable system ids\n";
  std::vector<unsigned int> stable_system_indices, unstable_system_indices;
  for(int i = 0; i < ens.nsys() ; i++ )
    {
      if(ens[i].is_disabled())
	unstable_system_indices.push_back(i);
      else
	stable_system_indices.push_back(i);
    }
  
  //  std::cerr << "# Writing stable system ids\n";
  if(cfg.count("stable_init_output"))
    {
      ofstream stable_init_output( cfg.optional("stable_init_output", string("stable_init_output.txt")).c_str() ); 
      print_selected_systems(ens_init,stable_system_indices, stable_init_output);
    }

  if(cfg.count("stable_final_output"))
    {
      ofstream stable_final_output( cfg.optional("stable_final_output", string("stable_final_output.txt")).c_str() ); 
      print_selected_systems(ens,stable_system_indices, stable_final_output);
    }

  if(cfg.count("unstable_init_output"))
    {
      ofstream unstable_init_output( cfg.optional("unstable_init_output", string("unstable_init_output.txt")).c_str() ); 
      print_selected_systems(ens_init,unstable_system_indices, unstable_init_output);
    }

  if(cfg.count("unstable_final_output"))
    {
      ofstream unstable_final_output( cfg.optional("unstable_final_output", string("unstable_final_output.txt")).c_str() ); 
      print_selected_systems(ens,unstable_system_indices, unstable_final_output);
    }

}

std::vector<std::vector<double> > calc_semimajor_axes(defaultEnsemble& ens)
{
  std::vector<std::vector<double> > semimajor_axes(ens.nsys(),std::vector<double>(ens.nbod(),0.));

  for(int sys_idx = 0; sys_idx < ens.nsys() ; sys_idx++)
    {
      int sys_id = ens[sys_idx].id();
      assert(sys_id>=0);
      assert(sys_id<ens.nsys());
      double star_mass = ens.mass(sys_idx,0);      
      double mass_effective = star_mass;
      double bx, by, bz, bvx, bvy, bvz;
      ens.get_barycenter(sys_idx,bx,by,bz,bvx,bvy,bvz,0);
      for(unsigned int bod=1;bod<ens.nbod();++bod) // Skip star since calculating orbits
	{
	  double mass = ens.mass(sys_idx,bod);
	  mass_effective += mass;
	  ens.get_barycenter(sys_idx,bx,by,bz,bvx,bvy,bvz,bod-1);
	  double x = ens.x(sys_idx,bod)-bx;
	  double y = ens.y(sys_idx,bod)-by;
	  double z = ens.z(sys_idx,bod)-bz;
	  double vx = ens.vx(sys_idx,bod)-bvx;
	  double vy = ens.vy(sys_idx,bod)-bvy;
	  double vz = ens.vz(sys_idx,bod)-bvz;
	  double a, e, i, O, w, M;
	  calc_keplerian_for_cartesian(a,e,i,O,w,M, x,y,z,vx,vy,vz, mass_effective);
	  semimajor_axes[sys_id][bod-1] = a;
	}
    }
  return semimajor_axes;
}

void disable_unstable_systems(defaultEnsemble& ens, const std::vector<std::vector<double> >& semimajor_axes_init, const double deltaa_threshold )
{
  for(int sys_idx = 0; sys_idx < ens.nsys() ; sys_idx++)
    {
      if(ens[sys_idx].is_disabled() ) continue;
      bool disable = false;
      int sys_id = ens[sys_idx].id();
      double star_mass = ens.mass(sys_idx,0);      
      double mass_effective = star_mass;
      double bx, by, bz, bvx, bvy, bvz;
      ens.get_barycenter(sys_idx,bx,by,bz,bvx,bvy,bvz,0);
      for(unsigned int bod=1;bod<ens.nbod();++bod) // Skip star since calculating orbits
	{
	  //	  std::cout << "body= " << bod << ": ";
	  //	  std::cout << "pos= (" << ens.x(sys_idx, bod) << ", " <<  ens.y(sys_idx, bod) << ", " << ens.z(sys_idx, bod) << ") vel= (" << ens.vx(sys_idx, bod) << ", " <<  ens.vy(sys_idx, bod) << ", " << ens.vz(sys_idx, bod) << ").\n";
	  double mass = ens.mass(sys_idx,bod);
	  mass_effective += mass;
	  ens.get_barycenter(sys_idx,bx,by,bz,bvx,bvy,bvz,bod-1);
	  double x = ens.x(sys_idx,bod)-bx;
	  double y = ens.y(sys_idx,bod)-by;
	  double z = ens.z(sys_idx,bod)-bz;
	  double vx = ens.vx(sys_idx,bod)-bvx;
	  double vy = ens.vy(sys_idx,bod)-bvy;
	  double vz = ens.vz(sys_idx,bod)-bvz;
	  double a, e, i, O, w, M;
	  calc_keplerian_for_cartesian(a,e,i,O,w,M, x,y,z,vx,vy,vz, mass_effective);
	  i *= 180/M_PI;
	  O *= 180/M_PI;
	  w *= 180/M_PI;
	  M *= 180/M_PI;

	  if(!((e>=0.)&&(e<1.0)))	{ disable = true; }
	  if(!((a>0.)&&(a<10.0)))	{ disable = true; }
	  assert(sys_id>=0);
	  assert(sys_id<semimajor_axes_init.size());
	  assert(bod-1>=0);
	  assert(bod-1<semimajor_axes_init[sys_id].size());
	  double da = a - semimajor_axes_init[sys_id][bod-1];
	  if(fabs(da) >= deltaa_threshold * semimajor_axes_init[sys_id][bod-1])
	    { disable = true; }

	  if(disable)
	    {
	      if(cfg.count("verbose"))
		std::cout << "# Disabling idx=" << sys_idx << " id=" << sys_id << " b=" << bod << " a= " << a << " e= " << e << " i= " << i << " Omega= " << O << " omega= " << w << " M= " << M << "\n";	  	
	      break;
	    }
	}
      if(disable) ens[sys_idx].set_disabled();
    }
}

bool needs_shrinking( const defaultEnsemble& ens ) 
{
  // This is the ratio we use when we are shrinking
  // if the ratio of active ones to all is less 
  // than this number then we trim the fat
  const double critical_ratio = 0.75;
  
  int count_disabled = number_of_disabled_systems( ens ) ;
  double ratio = double(count_disabled) / double( ens.nsys() );
  
  return ratio > critical_ratio;
}


/// Save a periodical snapshot if one is defined in the config file
/// Snapshot is usually saved as binary file, because it happens very
/// frequently and text files take very long time to generate.
void save_snapshot( defaultEnsemble& ens ) 
{
  if(cfg.count("snapshot")) snapshot::save( ens, cfg["snapshot"] );
}


/**
 *  This is a very crucial part of the Monte Carlo simulation
 *  We remove the disabled ones and make a smaller ensemble.
 *  We don't really need to make another ensemble. But keeping
 *  the same ensemble is a lot of trouble.
 *
 */
defaultEnsemble trim_disabled_systems( const defaultEnsemble& ens ) 
{
  int nsys = ens.nsys();
  int active_nsys = nsys - number_of_disabled_systems( ens ) ;
  // WARNING: Needed to add this to prevent segfaults.  TODO: Figure out why.
  if(active_nsys==0) return ens;  
  int nbod = ens.nbod();

  defaultEnsemble active_ens = defaultEnsemble::create( nbod, active_nsys );

  // Copy the active ones to the new ensemble
  for(int i = 0, j = 0; (i < nsys) && (j < active_nsys); i++)
    {
      if( !ens[i].is_disabled() )
	{
	  ens[i].copyTo( active_ens[j] );
	  j++;
	}
    }
  
  return active_ens;
}


void reactivate_systems(defaultEnsemble&ens)
{
  for(int i = 0; i < ens.nsys() ; i++)
    {
      if(ens[i].is_inactive())
	ens.set_active(i);
    }
}

volatile bool integration_loop_not_aborted_yet = true;
/**
 *   We can use this signal handler function
 *   to detect Ctrl-C and save the last snapshot and leave
 *   in a clean way
 *
 */
void ctrl_c_trap(int)
{
  fprintf(stderr, "# Break requested... Finishing the current GPU kernel call and will then save the results before exitting.\n");
  integration_loop_not_aborted_yet = false;
}
void catch_ctrl_c() 
{
  signal(SIGINT, &ctrl_c_trap );
}

int main(int argc, char* argv[] ) 
{
  // We keep it simple, later on one can use boost::program_options to 
  // have more options
  // but now we only use two configuration files. It is because the 
  // initial conditions configuration file can get really big and 
  // it has very little to do with other configuration options.
  if(argc < 3) cout << "Usage: montecarlo <integration configuration> <initial conditions configuration>" << endl;

  // First one is the configuration for integration
  string integ_configfile = argv[1];
  // the second one is the configuration for generating
  //  initial conditions it is used by \ref generate_ensemble_with_randomized_initial_conditions
  string initc_configfile = argv[2];
  
  cfg = config::load(integ_configfile);
  
  // 1.read keplerian coordinates from a file
  // 2.generate guesses based on the keplerian coordinates
  // 3.convert keplerian coordinates to an ensemble
  // The following line that is taken from swarm_tutorial_montecarlo.cpp
  // does the first three steps. Its pretty amazing.
  defaultEnsemble ens ; 
  if( cfg.count("input") ) 
    {    ens = snapshot::load(cfg["input"]);  }
  else
    {    ens = generate_ensemble_with_randomized_initial_conditions( config::load(initc_configfile) );  }
	
  // save the ensemble as a snapshot
  if(cfg.count("initial_snapshot"))
    {    snapshot::save( ens, cfg["initial_snapshot"] );	}

  defaultEnsemble ens_init = ens.clone() ;
  std::vector<std::vector<double> > semimajor_axes_init = calc_semimajor_axes(ens);

  // We want to find the ones that are really stable, so we integrate for
  // a really long time and over time we get rid of unstable ones. 
  double destination_time = cfg.optional("destination_time", 1.0E6);
  
  
  swarm::init(cfg);
  Pintegrator integ = integrator::create(cfg);
  integ->set_ensemble(ens);
  integ->set_destination_time ( destination_time );
  // We can set the following two items if we really need
  // longer integrations before we stop for checking the
  // ensemble and saving snapshots.
  // one kernel call to allow for prompt CPU pruning of unstable systems
  int max_kernel_calls_per_integrate = cfg.optional("max_kernel_calls_per_integrate",1);
  // 10^2 inner orbits at 200 time steps per inner orbit
  int max_itterations_per_kernel_call = cfg.optional("max_itterations_per_kernel_call",20000);
  integ->set_max_attempts( max_kernel_calls_per_integrate ); 
  integ->set_max_iterations ( max_itterations_per_kernel_call ); 
  SYNC;

  // integrate ensemble
  //  - drop the unstable ones as you go
  //  - redistribute the systems for better efficiency
  //  - save the results periodically
  //  This can be an infitie loop. but we can always get out of this
  //  the best way to do it is to use Ctrl-C. Further wecan
  //  define Ctrl-C signal handler to get out
  //  of this loop in a safe way. But we really want this loop
  //  to run for a long time
  reactivate_systems(ens);
  // EBF Experiment trying to expose host log.
  swarm::log::ensemble(*(swarm::log::manager::default_log()->get_hostlog()),ens);
  integ->flush_log();

  catch_ctrl_c();
  while( number_of_active_systems(ens) > 0 && integration_loop_not_aborted_yet ) {

    // 1. Integrate, we could use core_integrate but the general integrate
    // saves all the headache. We should only use core_integrate if we are
    // going to do everything on GPU, but since we are saving a snapshot in 
    // the middle, there's no point. It also has a nice for loop and can
    // to several kernel calls.
    integ->integrate();
    
    // 2. CPU-based tests to identify systems that can be terminated
    int active_ones = number_of_active_systems(ens);
    const double deltaa_frac_threshold = cfg.optional("deltaa_frac_threshold", 0.5);
    disable_unstable_systems( ens, semimajor_axes_init, deltaa_frac_threshold );
    std::cerr << "# Time: " << ens.time_ranges().min << " " << ens.time_ranges().median << " " << ens.time_ranges().max << " Systems: " << ens.nsys() << ", " << active_ones << ", ";
    active_ones = number_of_active_systems(ens);    
    std::cerr << active_ones << "\n";

    // EBF Experiment trying to expose host log.  
    swarm::log::ensemble_enabled(*(swarm::log::manager::default_log()->get_hostlog()),ens);
    integ->flush_log();

    // 3. Now we need to get rid of the inactive ones. There 
    // should be some criteria, whatever it is we are
    // going to put it in a function and call it \ref needs_shrinking
    if ( needs_shrinking( ens ) ) 
      {
	// Now this function will take care of trimming for us
	// We need to create a new ensemble of a smaller size
	// thats why we need to call set_ensemble again. because
	// the GPU ensemble should get recreated for us.
	//  we need better memory management to safely allow
	//  the following statement. Right now we don't have a
	//  very good automatic memory management
	//	std::cerr << "# Removing disabled systems\n";
	ens = trim_disabled_systems( ens );
	//	std::cerr << "# Testing if ens has any systems left.\n";
	if(ens.nsys()>0)
	  {
	    // std::cerr << "# Setting ensemble for integrator\n";
	    // WARNING: This appears to be the cause of SEGFAULT
	    //          if the ensemble is empty.  
	    integ->set_ensemble( ens );
	  }
	//	std::cerr << "# Time: " << ens.time_ranges() << " Total Systems: " << ens.nsys() << " Active Systems: " << active_ones << ", ";
	//	active_ones = number_of_active_systems(ens);    
	//	std::cerr << active_ones << "\n";
      }
    
    //    std::cerr << std::endl;

    // 4. We need to save a snapshot in case system crashes or someone
    // gets tired of it and hits Ctrl-C
    //    std::cerr << "# Saving snapshot\n";
    save_snapshot( ens );
    write_stable_systems(ens,ens_init);        
  }

  // Now we are at the end of the system, before we examine
  // the output we need to 
  //  std::cerr << "# Removing disabled systems\n";
  if(number_of_active_systems(ens)>0)
    {
      ens = trim_disabled_systems( ens );
      //  std::cerr << "# Saving snapshot\n";
      save_snapshot( ens );
    }

    write_stable_systems(ens,ens_init);    

}

