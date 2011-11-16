/*************
 *  Author : Saleh Dindar
 *
 *
 */
#include "common.hpp"
#include "utils.hpp"

using std::max;
using namespace swarm;
using std::string;

int number_of_disabled_systems(defaultEnsemble ens) {
	int count_running = 0;
	for(int i = 0; i < ens.nsys() ; i++)
		if( ens[i].is_disabled() ) count_running++;
	return count_running;
}

swarm::hostEnsemble generate_ensemble(swarm::config& cfg)  {
	int nsys = cfg.require("nsys",0);
	int nbod = cfg.require("nbod",0);
	double spacing_factor = cfg.optional( "spacing_factor", 1.4 );

	hostEnsemble ens = hostEnsemble::create( nbod, nsys );


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
			double rmag = pow( spacing_factor ,int(bod-1));  // semi-major axes exceeding this spacing results in systems are stable for nbody=3 and mass_planet=0.001
			double vmag = sqrt(mass_sun/rmag);  // spped for uniform circular motion
			double theta = (2.*M_PI*rand())/static_cast<double>(RAND_MAX);  // randomize initial positions along ecah orbit
			x  =  rmag*cos(theta); y  = rmag*sin(theta); z  = 0;
			vx = -vmag*sin(theta); vy = vmag*cos(theta); vz = 0.;

			// assign body a mass, position and velocity
			ens.set_body(sys, bod, mass_planet, x, y, z, vx, vy, vz);
		}
		ens[sys].set_active();
		ens[sys].time() = 0;
		ens[sys].id() = sys;
	}
	return ens;
}

double find_max_energy_conservation_error(ensemble& ens, ensemble& reference_ensemble ) {
	std::vector<double> energy_init(reference_ensemble.nsys());
	reference_ensemble.calc_total_energy(&energy_init[0]);
	std::vector<double> energy_final(ens.nsys());
	ens.calc_total_energy(&energy_final[0]);
	double max_deltaE = 0.;
	for(int sysid=0;sysid<ens.nsys();++sysid)
	{
	
		double deltaE = fabs ((energy_final[sysid]-energy_init[sysid])/energy_init[sysid] ) ;
		max_deltaE = max(deltaE, max_deltaE);
	}
	return max_deltaE;
}


bool validate_configuration(config& cfg){
  bool valid = true;                 // Indicates whether cfg parameters are valid
  int nsystems = cfg.optional("nsys",1000);
  int nbodypersystem = cfg.optional("nbod",3);
  // WARNING: blocksize isn't being used as the block size.  Trying to remove this.
//  int bs = cfg.optional("blocksize",16);

  // Check that parameters from command line are ok
//  if((bs<ENSEMBLE_CHUNK_SIZE)||(bs>64)) valid =false;
//  if( bs % ENSEMBLE_CHUNK_SIZE != 0 ) valid = false;
  if(!(nsystems>=1)||!(nsystems<=256000)) valid = false;
  if(!(nbodypersystem>=3)||!(nbodypersystem<=10)) valid = false;

  return valid;
}

void outputConfigSummary(std::ostream& o,swarm::config& cfg) {
	o << "# Integrator:\t" << cfg["integrator"] << "\n"
		<< "# Time step\t" << cfg["time_step"] << "\n"
		<< "# Destination time\t" << cfg["destination_time"] << "\n"
		<< "# Min time step\t" << cfg["min_time_step"] << "\n"
		<< "# Max time step\t" << cfg["max_time_step"] << "\n"
		<< "# No. Systems\t" << cfg["nsys"] << "\n"
		<< "# No. Bodies\t" << cfg["nbod"] << "\n"
//		<< "# Blocksize\t" << cfg["blocksize"] << "\n"
		<< std::endl;
}

// WARNING:  Is it really wise to provide all of these by default in utils?
//           This makes it easy to use a default value by accident
config default_config() {
	config cfg;
	cfg["nsys"] = 16;
	cfg["nbod"] = 3;
	cfg["integrator"] = "hermite"; // Set to use a GPU integrator
	cfg["time_step"] = "0.001";       // time step
	cfg["nbod"] = "3";
	cfg["nsys"] = "16";
//	cfg["blocksize"] = "16";
	cfg["log_writer"] = "null";
	return cfg;
}
std::ostream& operator << (std::ostream& o, const swarm::ensemble::range_t& r){
	if(r.min-r.max < 1e-11) 
		return o << r.median;
	else
		return o << r.median << "[" << (r.min-r.median) << "," << (r.max-r.median) << "] ";
}

bool compare_ensembles( swarm::ensemble& e1, swarm::ensemble &e2 , double & pos_diff, double & vel_diff, double & time_diff ) {
	if (e1.nsys() != e2.nsys() || e1.nbod() != e2.nbod() ) return false;

	pos_diff = vel_diff = time_diff = 0;

	for(int i = 0; i < e1.nsys(); i++) {
		for(int j = 0; j < e1.nbod() ; j++){

			double dp = sqrt( 
					  square ( e1[i][j][0].pos() - e2[i][j][0].pos() ) 
					+ square ( e1[i][j][1].pos() - e2[i][j][1].pos() ) 
					+ square ( e1[i][j][2].pos() - e2[i][j][2].pos() ) ) ;

			double dv = sqrt( 
					  square ( e1[i][j][0].vel() - e2[i][j][0].vel() ) 
					+ square ( e1[i][j][1].vel() - e2[i][j][1].vel() ) 
					+ square ( e1[i][j][2].vel() - e2[i][j][2].vel() ) ) ;

			if ( dp > pos_diff ) pos_diff = dp;
			if ( dv > vel_diff ) vel_diff = dv;

		}

		double dt = fabs(e1[i].time() - e2[i].time());
		if ( dt > time_diff ) time_diff = dt;

	}
	return true;
}

