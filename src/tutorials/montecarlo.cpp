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
#include "random.hpp"
#include "kepler.hpp"

#define SYNC cudaThreadSynchronize()


using namespace swarm;
using namespace std;

config cfg;

void inspect(defaultEnsemble &ens, const int sys, const int bod ) {
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
defaultEnsemble generate_randomized_initial_conditions(const config& cfg) 
{
	defaultEnsemble ens = defaultEnsemble::create( cfg.require("nbod",0), cfg.require("nsys",0) );

  std::cerr << "Set initial time for all systems = ";
  double time_init = cfg.optional("time_init", 0.0);
  std::cerr << time_init << ".\n";

  bool use_jacobi = cfg.optional("use_jacobi", 0);

  for(unsigned int sys=0;sys<ens.nsys();++sys)
    {
		ens[sys].id() = sys;
		ens[sys].time() = time_init;
		ens[sys].set_active();
      // set sun to unit mass and at origin
      double mass_star = cfg.optional("mass_star", 1.);
      double x=0, y=0, z=0, vx=0, vy=0, vz=0;
      ens.set_body(sys, 0, mass_star, x, y, z, vx, vy, vz);

      double mass_enclosed = mass_star;
      for(unsigned int bod=1;bod<ens.nbod();++bod)
	{
	  double mass_planet = draw_value_from_config(cfg,"mass",bod,0.,mass_star);
	  mass_enclosed += mass_planet;

	  double a = draw_value_from_config(cfg,"a",bod,0.001,10000.);
	  double e = draw_value_from_config(cfg,"ecc",bod,0.,1.);
	  double i = draw_value_from_config(cfg,"inc",bod,-180.,180.);
	  double O = draw_value_from_config(cfg,"node",bod,-720.,720.);
	  double w = draw_value_from_config(cfg,"omega",bod,-720.,720.);
	  double M = draw_value_from_config(cfg,"meananom",bod,-720.,720.);
//	  	  std::cout << "# Drawing sys= " << sys << " bod= " << bod << ' ' << mass_planet << "  " << a << ' ' << e << ' ' << i << ' ' << O << ' ' << w << ' ' << M << '\n';

	  i *= M_PI/180.;
	  O *= M_PI/180.;
	  w *= M_PI/180.;
	  M *= M_PI/180.;

	  double mass = use_jacobi ? mass_enclosed : mass_star+mass_planet;
	  calc_cartesian_for_ellipse(x,y,z,vx,vy,vz, a, e, i, O, w, M, mass);
	   
	//  printf("%d %d: %lg (%lg %lg %lg) (%lg %lg %lg) \n", sys, bod, mass_planet, x,y,z,vx,vy,vz);

	  double bx, by, bz, bvx, bvy, bvz;
	  ens.get_barycenter(sys,bx,by,bz,bvx,bvy,bvz,bod-1);
	  x  += bx;	  y  += by;	  z  += bz;
	  vx += bvx;	  vy += bvy;	  vz += bvz;

	  // assign body a mass, position and velocity
	  ens.set_body(sys, bod, mass_planet, x, y, z, vx, vy, vz);
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
  return ens;
}


/**
 * This function uses DEFINES. Now I personally hate defines. Because sometimes
 * you wanna do things runtime and then it is not possible. Why do you have to use
 * define. Just use a normal if statement and a global variable that tells you
 * which one to use. what if I define JACOBI and BARRYCENTRIC at the same time.
 * And the way that they misspelled BARRYCENTRIC. Sounds like all planets revolve
 * around a guy named BARRY.
 *
 */
void print_selected_systems_for_demo(swarm::ensemble& ens)
{
	enum {
		JACOBI, BARYCENTRIC, ASTROCENTRIC
	} COORDINATE_SYSTEM = BARYCENTRIC;


  std::streamsize cout_precision_old = std::cout.precision();
  std::cout.precision(10);
  unsigned int nprint = std::min(10,ens.nsys());
  for(unsigned int systemid = 0; systemid< nprint; ++systemid)
    {
      std::cout << "sys= " << systemid << " time= " << ens.time(systemid) << "\n";
      double star_mass = ens.mass(systemid,0);
      double mass_effective = star_mass;
      double bx, by, bz, bvx, bvy, bvz;
	  switch(COORDINATE_SYSTEM) {
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
	  std::cout << "body= " << bod << ": ";
	  //	  std::cout << "pos= (" << ens.x(systemid, bod) << ", " <<  ens.y(systemid, bod) << ", " << ens.z(systemid, bod) << ") vel= (" << ens.vx(systemid, bod) << ", " <<  ens.vy(systemid, bod) << ", " << ens.vz(systemid, bod) << ").\n";
	  double mass = ens.mass(systemid,bod);
	  	
	  	switch(COORDINATE_SYSTEM) {
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
	      std::cout << " a= " << a << " e= " << e << " i= " << i << " Omega= " << O << " omega= " << w << " M= " << M << "\n";


	}
    }
  std::cout.precision(cout_precision_old);
  std::cout << std::flush;
}

bool needs_shrinking( const defaultEnsemble& ens ) {
	// This is the ratio we use when we are shrinking
	// if the ratio of active ones to all is less 
	// than this number then we trim the fat
	const double critical_ratio = .5;

	int count_disabled = number_of_disabled_systems( ens ) ;
	double ratio = double(count_disabled) / double( ens.nsys() );

	return ratio > critical_ratio;
}


/// Save a periodical snapshot if one is defined in the config file
/// Snapshot is usually saved as binary file, because it happens very
/// frequently and text files take very long time to generate.
void save_snapshot( defaultEnsemble& ens ) {
	if(cfg.count("snapshot")) snapshot::save( ens, cfg["snapshot"] );
}


/**
 *  This is a very crucial part of the Monte Carlo simulation
 *  We remove the disabled ones and make a smaller ensemble.
 *  We don't really need to make another ensemble. But keeping
 *  the same ensemble is a lot of trouble.
 *
 */
defaultEnsemble trim_disabled_systems( const defaultEnsemble& ens ) {
	int nsys = ens.nsys();
	int active_nsys = nsys - number_of_disabled_systems( ens ) ;
	int nbod = ens.nbod();

	defaultEnsemble active_ens = defaultEnsemble::create( nbod, active_nsys );

	// Copy the active ones to the new ensemble
	for(int i = 0, j = 0; i < nsys && j < active_nsys; i++)
		if( !ens[i].is_disabled() ) 
			ens[i].copyTo( active_ens[j] ), j++;

	return active_ens;
}

void reactivate_systems(defaultEnsemble&ens){
		for(int i = 0; i < ens.nsys() ; i++)
			if(ens[i].state() == ensemble::Sys::SYSTEM_INACTIVE )
				ens[i].set_active();
}
volatile bool integration_loop_not_aborted_yet = true;
/**
 *   We can use this signal handler function
 *   to detect Ctrl-C and save the last snapshot and leave
 *   in a clean way
 *
 */
void ctrl_c_trap(int){
	fprintf(stderr, "Break requested, saving the results\n");
	integration_loop_not_aborted_yet = false;
}
void catch_ctrl_c() {
	signal(SIGINT, &ctrl_c_trap );
}

int main(int argc, char* argv[] ) {
	// We keep it simple, later on one can use boost::program_options to 
	// have more options
	// but now we only use two configuration files. It is because the 
	// initial conditions configuration file can get really big and 
	// it has very little to do with other configuration options.
	if(argc < 3) cout << "Usage: montecarlo <integration configuration> <initial conditions configuration>" << endl;

	// First one is the configuration for integration
	string integ_configfile = argv[1];
	// the second one is the configuration for generating
	//  initial conditions it is used by \ref generate_randomized_initial_conditions
	string initc_configfile = argv[2];

	cfg = config::load(integ_configfile);

	// 1.read keplerian coordinates from a file
	// 2.generate guesses based on the keplerian coordinates
	// 3.convert keplerian coordinates to an ensemble
	// The following line that is taken from swarm_tutorial_montecarlo.cpp
	// does the first three steps. Its pretty amazing.
	defaultEnsemble ens ; 
	if( cfg.count("input") ) {
		ens = snapshot::load(cfg["input"]);
	}else{
		ens = generate_randomized_initial_conditions( config::load(initc_configfile) );
	}
	
	// save the ensemble as a snapshot
	if(cfg.count("initial_snapshot")){
		snapshot::save( ens, cfg["initial_snapshot"] );	
	}

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
	integ->set_max_attempts( 5 );
	// integ->set_max_iterations ( ? );
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
	catch_ctrl_c();
	while( number_of_active_systems(ens) > 0 && integration_loop_not_aborted_yet ) {

		// 1. Integrate, we could use core_integrate but the general integrate
		// saves all the headache. We should only use core_integrate if we are
		// going to do everything on GPU, but since we are saving a snapshot in 
		// the middle, there's no point. It also has a nice for loop and can
		// to several kernel calls.
		integ->integrate();

		int active_ones = number_of_active_systems(ens);

		cout << ens.time_ranges() << ", " << active_ones << endl;


		// 2. Now we need to get rid of the inactive ones. There 
		// should be some criteria, whatever it is we are
		// going to put it in a function and call it \ref needs_shrinking
		if ( needs_shrinking( ens ) ) {
			// Now this function will take care of trimming for us
			// We need to create a new ensemble of a smaller size
			// thats why we need to call set_ensemble again. because
			// the GPU ensemble should get recreated for us.
			//  we need better memory management to safely allow
			//  the following statement. Right now we don't have a
			//  very good automatic memory management
			ens = trim_disabled_systems( ens );
			integ->set_ensemble( ens );
		}

		// 3. We need to save a snapshot in case system crashes or someone
		// gets tired of it and hits Ctrl-C
		save_snapshot( ens );

	}

	// Now we are at the end of the system, before we examine
	// the output we need to 
	ens = trim_disabled_systems( ens );
	save_snapshot( ens );
	
	// find the stable ones and output the initial conditions for the stable
	// ones in keplerian coordinates
	ofstream output( cfg.optional("output", string("output.txt")).c_str() ); 

	for(int i = 0; i < ens.nsys() ; i++ ){
		if(!ens[i].is_disabled())
			output << ens[i].id() << endl;
	}
	
}

