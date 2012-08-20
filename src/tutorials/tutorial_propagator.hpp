/*
 *  This is a simple tutorial used in doxygen pages
 *  should go through program2doxygen before it
 *  can be used by doxygen.
 *  
 *
 */
// \page TutorialPropagator Tutorial for making a Propagator
//
// This is the header file where you define the new propagator
//
#include "swarm/swarmplugin.h"


struct TutorialPropagatorParams {
	double time_step;
	TutorialPropagatorParams(const config& cfg){
		time_step = cfg.require("time_step", 0.0);
	}
};

template<class T,class Gravitation>
struct TutorialPropagator {
	typedef TutorialPropagatorParams params;
	const static int nbod = T::n;

	params _params;


	// Runtime variables
	ensemble::SystemRef& sys;
	Gravitation& calcForces;
	int b;
	int c;
	int ij;
	bool body_component_grid;
	bool first_thread_in_system;
	double max_timestep;


	GPUAPI TutorialPropagator(const params& p,ensemble::SystemRef& s,
			Gravitation& calc)
		:_params(p),sys(s),calcForces(calc){}

	GPUAPI void init()  { }

	GPUAPI void shutdown() { }

        GPUAPI void convert_internal_to_std_coord() {} 
        GPUAPI void convert_std_to_internal_coord() {}

	static GENERIC int thread_per_system(){
		return nbod * 3;
	}

	static GENERIC int shmem_per_system() {
		 return 0;
	}

	GPUAPI bool is_in_body_component_grid()
        { return  ((b < nbod) && (c < 3)); }	

	GPUAPI bool is_in_body_component_grid_no_star()
        { return ( (b!=0) && (b < nbod) && (c < 3) ); }	

	GPUAPI bool is_first_thread_in_system()
        { return (thread_in_system()==0); }	

	GPUAPI void advance(){
		double h = min(_params.time_step, max_timestep);
		double pos = 0.0, vel = 0.0;
		double acc = 0.0, jerk = 0.0;
		const double third = 1.0/3.0;
		
		if( is_in_body_component_grid() )
			pos = sys[b][c].pos() , vel = sys[b][c].vel();


		calcForces(ij,b,c,pos,vel,acc,jerk);
		// Integrator
		pos = pos +  h*(vel+(h*0.5)*(acc+(h*third)*jerk));
		vel = vel +  h*(acc+(h*0.5)*jerk);


		// Finalize the step
		if( is_in_body_component_grid() )
			sys[b][c].pos() = pos , sys[b][c].vel() = vel;
		if( is_first_thread_in_system() ) 
			sys.time() += h;
	}
};
