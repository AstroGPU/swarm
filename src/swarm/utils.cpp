/*************************************************************************
 * Copyright (C) 2011 by Saleh Dindar and the Swarm-NG Development Team  *
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

/*! \file util.cpp
 *    \brief Implements the utility functions for swarm. 
 *
 * @authors Saleh Dindar
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


/**
 *
 * The velocity vector for planets is perpendicular to the position vector
 * (from star to planet) and the magnitude depends on the ejection_factor
 * |v| = sqrt(2*G*M/r)*ejection_factor. 
 * We get different type of orbit based on ejection_factor:
 *    1/sqrt(2) : circular orbit
 *    < 1 : elliptical orbit
 *    = 1 : parabolic orbit
 *    > 1 : hyperbolic orbit
 *
 * Configuration options:
 *   nsys: Number of systems in the ensemble
 *   nbod: Number of bodies per system
 *   spacing_factor: determines the spacing between
 *     the planets, distance of planet i from star is spacing_factor 
 *     times the distance of planet i-1 from star.
 *   ejection_factor: determines the type of orbit see above
 *   planet_mass: ratio of the planet mass to the star mass. defaults
 *   to Jupiter mass planets (0.001)
 *
 *
 */
swarm::hostEnsemble generate_ensemble(swarm::config& cfg)  {
	int nsys = cfg.require("nsys",0);
	int nbod = cfg.require("nbod",0);
	double spacing_factor = cfg.optional( "spacing_factor", 1.4 );
    double planet_mass = cfg.optional( "planet_mass" , .001 );
	double ejection_factor  = cfg.optional("ejection_factor", 1.0/sqrt(2) );


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
			double rmag = pow( spacing_factor ,int(bod-1));  // semi-major axes exceeding this spacing results in systems are stable for nbody=3 and mass_planet=0.001
			double vmag = sqrt(2*mass_sun/rmag);  // spped for uniform circular motion
			double theta = (2.*M_PI*rand())/static_cast<double>(RAND_MAX);  // randomize initial positions along ecah orbit
			x  =  rmag*cos(theta); y  = rmag*sin(theta); z  = 0;
			vx = -vmag*sin(theta); vy = vmag*cos(theta); vz = 0.;

			// assign body a mass, position and velocity
			ens.set_body(sys, bod, planet_mass , x, y, z, vx, vy, vz);
		}
		ens[sys].set_active();
		ens[sys].time() = 0;
		ens[sys].id() = sys;
	}
	return ens;
}

double find_max_energy_conservation_error(ensemble& ens, ensemble& reference_ensemble ) {
    return energy_conservation_error_range(ens,reference_ensemble).max;
}
ensemble::range_t energy_conservation_error_range(ensemble& ens, ensemble& reference_ensemble ) {
	std::vector<double> 
            energy_init(reference_ensemble.nsys())
            ,energy_final(ens.nsys())
            ,deltaE(ens.nsys());

	reference_ensemble.calc_total_energy(&energy_init[0]);
	ens.calc_total_energy(&energy_final[0]);

	for(int sysid=0;sysid<ens.nsys();++sysid)
		deltaE[sysid] =  fabs ((energy_final[sysid]-energy_init[sysid])/energy_init[sysid] ) ;

        return ensemble::range_t::calculate( deltaE.begin(), deltaE.end() );
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
  if(!(nbodypersystem>=3)||!(nbodypersystem<=MAX_NBODIES)) valid = false;

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

			// Distance between body position in ensemble e1 and e2
			double dp = sqrt( 
					  square ( e1[i][j][0].pos() - e2[i][j][0].pos() ) 
					+ square ( e1[i][j][1].pos() - e2[i][j][1].pos() ) 
					+ square ( e1[i][j][2].pos() - e2[i][j][2].pos() ) ) ;

			// Average magnitude of the position in e1 and e2
			double ap = sqrt( 
					  square ( e1[i][j][0].pos() + e2[i][j][0].pos() ) 
					+ square ( e1[i][j][1].pos() + e2[i][j][1].pos() ) 
					+ square ( e1[i][j][2].pos() + e2[i][j][2].pos() ) ) / 2.0 ;

			// Difference between body velocities in ensemble e1 and e2
			double dv = sqrt( 
					  square ( e1[i][j][0].vel() - e2[i][j][0].vel() ) 
					+ square ( e1[i][j][1].vel() - e2[i][j][1].vel() ) 
					+ square ( e1[i][j][2].vel() - e2[i][j][2].vel() ) ) ;

			// Average magnitude of the velocity in e1 and e2
			double av = sqrt( 
					  square ( e1[i][j][0].vel() + e2[i][j][0].vel() ) 
					+ square ( e1[i][j][1].vel() + e2[i][j][1].vel() ) 
					+ square ( e1[i][j][2].vel() + e2[i][j][2].vel() ) ) / 2.0 ;

			if ( dp > pos_diff ) pos_diff = dp;
			if ( dv > vel_diff ) vel_diff = dv;

		}

		// Difference between time divided by average of the two
		double dt = fabs(e1[i].time() - e2[i].time())
			/(e1[i].time() + e2[i].time())*2;
		if ( dt > time_diff ) time_diff = dt;

	}
	return true;
}

