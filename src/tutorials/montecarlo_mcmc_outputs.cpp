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

/**
 * Read in Keplerian coordinates from a text file
 *
 */
defaultEnsemble generate_ensemble_with_initial_conditions_keplerian_from_file(const config& cfg) 
{
  defaultEnsemble ens = defaultEnsemble::create( cfg.require("nbod",0), cfg.require("nsys",0) );
  std::cerr << "# nsystems = " << ens.nsys() << " nbodies/system = " << ens.nbod() << "\n";

  //  std::cerr << "# Set initial time for all systems = ";
  double time_init = cfg.optional("time_init", 0.0);
  //  std::cerr << time_init << ".\n";

  const bool use_jacobi = cfg.optional("use_jacobi", 0);

  //  std::cerr << "# Writing stable system ids\n";
  ifstream jacobi_input( cfg.optional("input_mcmc_keplerian", string("mcmc.out")).c_str() ); 

  for(unsigned int sysid=0;sysid<ens.nsys();++sysid)
    {
      assert(jacobi_input.good());
      string id;
      double mass_star;
      jacobi_input >> id >> mass_star;
      ens[sysid].id() = sysid;
      ens[sysid].time() = time_init;
      ens[sysid].set_active();
      double x=0, y=0, z=0, vx=0, vy=0, vz=0;
      ens.set_body(sysid, 0, mass_star, x, y, z, vx, vy, vz);
      double mass_enclosed = mass_star;
      for(unsigned int bod=1;bod<ens.nbod();++bod)
	{

	  double mass_planet, a, e, i, O, w, M;
	  jacobi_input >> mass_planet >> a >> e >> i >> w >> O >> M;
	  //	  a = pow((period/365.25)*(period/365.25)*mass_enclosed,1.0/3.0);
	  i *= M_PI/180.;
	  O *= M_PI/180.;
	  w *= M_PI/180.;
	  M *= M_PI/180.;

	  mass_enclosed += mass_planet;	  
	  double mass = use_jacobi ? mass_enclosed : mass_star+mass_planet;
	  if(cfg.count("verbose")&&(sysid<10))
	    std::cout << "# Drawing sysid= " << sysid << " bod= " << bod << ' ' << mass_planet << "  " << a << ' ' << e << ' ' << i*180./M_PI << ' ' << O*180./M_PI << ' ' << w*180./M_PI << ' ' << M*180./M_PI << '\n';
	  
	  calc_cartesian_for_ellipse(x,y,z,vx,vy,vz, a, e, i, O, w, M, mass);
      //  printf("%d %d: %lg (%lg %lg %lg) (%lg %lg %lg) \n", sysid, bod, mass_planet, x,y,z,vx,vy,vz);

	  if(use_jacobi)
	    {
	      double bx, by, bz, bvx, bvy, bvz;
	      ens.get_barycenter(sysid,bx,by,bz,bvx,bvy,bvz,bod-1);
	      x  += bx;	  y  += by;	  z  += bz;
	      vx += bvx;  vy += bvy;	  vz += bvz;
	    }
	  
	  // assign body a mass, position and velocity
	  ens.set_body(sysid, bod, mass_planet, x, y, z, vx, vy, vz);

	  if(cfg.count("verbose")&&(sysid<10))
	    {
	      double x_t = ens.x(sysid,bod);
	      double y_t = ens.y(sysid,bod);
	      double z_t = ens.z(sysid,bod);
	      double vx_t = ens.vx(sysid,bod);
	      double vy_t = ens.vy(sysid,bod);
	      double vz_t = ens.vz(sysid,bod);
	      
	      std::cout << " x= " << x << "=" << x_t << " ";
	      std::cout << " y= " << y << "=" << y_t << " ";
	      std::cout << " z= " << z << "=" << z_t << " ";
	      std::cout << "vx= " << vx << "=" << vx_t << " ";
	      std::cout << "vy= " << vy << "=" << vy_t << " ";
	      std::cout << "vz= " << vz << "=" << vz_t << "\n";
	    }
      
	}  // end loop over bodies
  
      // Shift into barycentric frame
      ens.get_barycenter(sysid,x,y,z,vx,vy,vz);
      for(unsigned int bod=0;bod<ens.nbod();++bod)
	{
	  ens.set_body(sysid, bod, ens.mass(sysid,bod), 
		       ens.x(sysid,bod)-x, ens.y(sysid,bod)-y, ens.z(sysid,bod)-z, 
		       ens.vx(sysid,bod)-vx, ens.vy(sysid,bod)-vy, ens.vz(sysid,bod)-vz);	  
	}  // end loop over bodies

    } // end loop over systems
  return ens;
}


/**
 * Read in Cartesian coordinates from a text file
 *
 */
defaultEnsemble generate_ensemble_with_initial_conditions_cartesian_from_file(const config& cfg) 
{
  defaultEnsemble ens = defaultEnsemble::create( cfg.require("nbod",0), cfg.require("nsys",0) );
  std::cerr << "# nsystems = " << ens.nsys() << " nbodies/system = " << ens.nbod() << "\n";

  //  std::cerr << "# Set initial time for all systems = ";
  double time_init = cfg.optional("time_init", 0.0);
  //  std::cerr << time_init << ".\n";

  const bool use_jacobi = cfg.optional("use_jacobi", 0);

  //  std::cerr << "# Writing stable system ids\n";
  ifstream mcmc_input( cfg.optional("input_mcmc_cartesian", string("mcmc.out")).c_str() ); 

  assert(mcmc_input.good());

  const int lines_to_skip = cfg.optional("lines_to_skip", 0);
  for(int i=0;i<lines_to_skip;++i)
    {
      assert(mcmc_input.good());
      char junk[1024];
      mcmc_input.getline(junk,1024);
    }

  for(unsigned int sysid=0;sysid<ens.nsys();++sysid)
    {
      assert(mcmc_input.good());
#if 1
      std::vector<double> masses(ens.nbod(),0.0);
      for(unsigned int bod=0;bod<ens.nbod();++bod)
	{
	  mcmc_input >> masses[bod];
	}  // end loop over bodies

      ens[sysid].id() = sysid;
      ens[sysid].time() = time_init;
      ens[sysid].set_active();
      double x=0, y=0, z=0, vx=0, vy=0, vz=0;
      ens.set_body(sysid, 0, masses[0], x, y, z, vx, vy, vz);
      double mass_enclosed = masses[0];
      for(unsigned int bod=1;bod<ens.nbod();++bod)
	{
	  mcmc_input >> x >> y >> z >> vx >> vy >> vz;

	  mass_enclosed += masses[bod];
	  double mass = use_jacobi ? mass_enclosed : masses[0]+masses[bod];
	  if(cfg.count("verbose")&&(sysid<10))
	    std::cout << "# Drawing sysid= " << sysid << " bod= " << bod << ' ' << masses[bod] << "  " << x << ' ' << y << ' ' << z << ' ' << vx << ' ' << vy << ' ' << vz << '\n';
	  
	  if(use_jacobi)
	    {
	      double bx, by, bz, bvx, bvy, bvz;
	      ens.get_barycenter(sysid,bx,by,bz,bvx,bvy,bvz,bod-1);
	      x  += bx;	  y  += by;	  z  += bz;
	      vx += bvx;  vy += bvy;	  vz += bvz;
	    }
	  
	  // assign body a mass, position and velocity
	  ens.set_body(sysid, bod, masses[bod], x, y, z, vx, vy, vz);

	  if(cfg.count("verbose")&&(sysid<10))
	    {
	      double x_t = ens.x(sysid,bod);
	      double y_t = ens.y(sysid,bod);
	      double z_t = ens.z(sysid,bod);
	      double vx_t = ens.vx(sysid,bod);
	      double vy_t = ens.vy(sysid,bod);
	      double vz_t = ens.vz(sysid,bod);
	      
	      std::cout << " x= " << x << "=" << x_t << " ";
	      std::cout << " y= " << y << "=" << y_t << " ";
	      std::cout << " z= " << z << "=" << z_t << " ";
	      std::cout << "vx= " << vx << "=" << vx_t << " ";
	      std::cout << "vy= " << vy << "=" << vy_t << " ";
	      std::cout << "vz= " << vz << "=" << vz_t << "\n";
	    }
#else
      std::vector<double> masses(ens.nbod(),0.0);
      double x=0, y=0, z=0, vx=0, vy=0, vz=0;
      double mass_enclosed = 0.0;
      ens[sysid].id() = sysid;
      ens[sysid].time() = time_init;
      ens[sysid].set_active();
      for(unsigned int bod=0;bod<ens.nbod();++bod)
	{
	  mcmc_input >> masses[bod];
	  mcmc_input >> x >> y >> z >> vx >> vy >> vz;

	  mass_enclosed += masses[bod];
	  double mass = use_jacobi ? mass_enclosed : masses[0]+masses[bod];
	  if(cfg.count("verbose")&&(sysid<10))
	    std::cout << "# Drawing sysid= " << sysid << " bod= " << bod << ' ' << masses[bod] << "  " << x << ' ' << y << ' ' << z << ' ' << vx << ' ' << vy << ' ' << vz << '\n';
	  
	  if(use_jacobi)
	    {
	      double bx, by, bz, bvx, bvy, bvz;
	      ens.get_barycenter(sysid,bx,by,bz,bvx,bvy,bvz,bod-1);
	      x  += bx;	  y  += by;	  z  += bz;
	      vx += bvx;  vy += bvy;	  vz += bvz;
	    }
	  
	  // assign body a mass, position and velocity
	  ens.set_body(sysid, bod, masses[bod], x, y, z, vx, vy, vz);

	  if(cfg.count("verbose")&&(sysid<10))
	    {
	      double x_t = ens.x(sysid,bod);
	      double y_t = ens.y(sysid,bod);
	      double z_t = ens.z(sysid,bod);
	      double vx_t = ens.vx(sysid,bod);
	      double vy_t = ens.vy(sysid,bod);
	      double vz_t = ens.vz(sysid,bod);
	      
	      std::cout << " x= " << x << "=" << x_t << " ";
	      std::cout << " y= " << y << "=" << y_t << " ";
	      std::cout << " z= " << z << "=" << z_t << " ";
	      std::cout << "vx= " << vx << "=" << vx_t << " ";
	      std::cout << "vy= " << vy << "=" << vy_t << " ";
	      std::cout << "vz= " << vz << "=" << vz_t << "\n";
	    }
#endif      
=======
      
>>>>>>> Stashed changes
	}  // end loop over bodies
  
      // Shift into barycentric frame
      ens.get_barycenter(sysid,x,y,z,vx,vy,vz);
      for(unsigned int bod=0;bod<ens.nbod();++bod)
	{
	  ens.set_body(sysid, bod, ens.mass(sysid,bod), 
		       ens.x(sysid,bod)-x, ens.y(sysid,bod)-y, ens.z(sysid,bod)-z, 
		       ens.vx(sysid,bod)-vx, ens.vy(sysid,bod)-vy, ens.vz(sysid,bod)-vz);	  
	}  // end loop over bodies

    } // end loop over systems
  return ens;
}


void print_system(const swarm::ensemble& ens, const int systemid, std::ostream &os = std::cout)
{
  enum {
    JACOBI, BARYCENTRIC, ASTROCENTRIC
  } COORDINATE_SYSTEM = BARYCENTRIC;
  const bool use_jacobi = cfg.optional("use_jacobi_output", 0);  
  if(use_jacobi_output) COORDINATE_SYSTEM = JACOBI;

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
		std::cout << "# Disabling idx=" << sys_idx << " id=" << sys_id << " b=" << bod << " ainit= " << semimajor_axes_init[sys_id][bod-1] << " a= " << a << " e= " << e << " i= " << i << " Omega= " << O << " omega= " << w << " M= " << M << "\n";	  	
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
  if(argc < 3) cout << "Usage: montecarlo_mcmc_outputs <integration configuration> <initial conditions configuration>" << endl;

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
  else if(cfg.count("input_mcmc_keplerian") )
    {    ens = generate_ensemble_with_initial_conditions_keplerian_from_file( config::load(initc_configfile) );  }
  else if(cfg.count("input_mcmc_cartesian") )
    {    ens = generate_ensemble_with_initial_conditions_cartesian_from_file( config::load(initc_configfile) );  }
  else
    {
      std::cerr << "# Must specify one of input [for binary snapshot], input_mcmc_keplerian or input_mcmc_cartesian.\n";
      return 255;
    }
	
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
  int num_integrate_calls = 0;
  while( number_of_active_systems(ens) > 0 && integration_loop_not_aborted_yet ) {

    // 1. Integrate, we could use core_integrate but the general integrate
    // saves all the headache. We should only use core_integrate if we are
    // going to do everything on GPU, but since we are saving a snapshot in 
    // the middle, there's no point. It also has a nice for loop and can
    // to several kernel calls.
    integ->integrate();
    ++num_integrate_calls;

    // 2. CPU-based tests to identify systems that can be terminated
    int active_ones = number_of_active_systems(ens);
    const double deltaa_frac_threshold = cfg.optional("deltaa_frac_threshold", 0.5);
    disable_unstable_systems( ens, semimajor_axes_init, deltaa_frac_threshold );
    double max_deltaE = find_max_energy_conservation_error(ens, ens_init );
    std::cerr << "# Time: " << ens.time_ranges().min << " " << ens.time_ranges().median << " " << ens.time_ranges().max << " Systems: " << ens.nsys() << ", " << active_ones << ", ";
    active_ones = number_of_active_systems(ens);    
    std::cerr << active_ones << "  max|dE/E|= " << max_deltaE << "\n";

    // EBF Experiment trying to expose host log.  
    if(num_integrate_calls%10==0)
      {
	swarm::log::ensemble_enabled(*(swarm::log::manager::default_log()->get_hostlog()),ens);
	integ->flush_log();
      }

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

