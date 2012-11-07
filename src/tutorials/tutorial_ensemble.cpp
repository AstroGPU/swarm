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

/*! \file tutorial_ensemble.cpp
 *   \brief A tutorial on how to manipulate ensembles
 * 
 * 
 */
// @page TutorialEnsemble Tutorial on using the Ensemble data structure.
// In this tutorial, we demonstrate how to generate an ensemble from
// a mathematical description.
//
// We will create an ensemble object, fill it with information of
// planetary systems that we synthesize. In the end we save the whole
// ensemble to the file. However, one can go ahead and do a number of
// arbitrary integrations on the ensemble object as demonstrated in
// @ref TutorialBeginner.
//
// We generate planetary systems where all the planets are on a parabola
// The star is placed at the origin. and the planets are on 
// y = x^2/4/R - R. 
//
// Let's get into the code.
// First some definitions and includes.
#include "swarm/swarm.h"
#include <iostream>
using namespace swarm;
using namespace std;

// Here we define some properties to be constant, like :
// 
// Number of bodies in a system
const int nbod = 3;
// Number of planetary systems
const int nsys = 64;
// Mass of the stars and the planets
const double mass_star   = 1;
const double mass_planet = 0.0001;
// X component of the planet positions.
const double xpos[4] = {  -2.0, +2.0, -4.0, +4.0 };
// mu = GM , but we have set both G and M to 1.
const double mu = 1; 
// the minimum distance of the planets to the star
// on the parabolic orbit. or in other terms, the distance between
// the bottom of the parabola and the star (at origin).
const double R = 1; 
//
// We create a convenience norm function
double norm(double x,double y){
	return sqrt(x*x+y*y);
}
//
// Now we can start the main procedure.
// 
int main(int argc, char* argv[]){
//  our program take one ARGV parameter: the
// name of the output file.
if(argc <= 1){
	cout << "Usage: parabolic_collision <outputfilename>" << endl;
}
const string outputfn = argv[1];

// First we have to initialize swarm, we don't need any specific configurations
// since we are not performing any integrations.
init(config());
// Create an empty ensemble.
defaultEnsemble ens = defaultEnsemble::create(nbod,nsys);

// Now we have to populate the ensemble with data, we have
// to write two nesting for loops, one over systems
// and one over planets.
for(int i = 0; i < nsys ; i++){
	ensemble::SystemRef s = ens[i];
	
	// Set the system parameters, the system must be set active, otherwise
	// it will not be integrated.
	s.id() = 0;
	s.time() = 0;
	s.set_active();
	
	// Place a stationary star at origin
	s[0].mass() = mass_star;
	s[0].x()  = 0, s[0].y()  = 0, s[0].z()  = 0 ;
	s[0].vx() = 0, s[0].vy() = 0, s[0].vz() = 0 ;
	
	
	// Place the rest of the planets on the parabola
	for(int b = 1; b < nbod; b++){
		// x comes from the array (defined at the top). Y is calculated
		// to be on the parabola.
		double x = xpos[b-1], y = x*x/4/R - R ;
		// Speed is calculated based on orbital parameters for a parabolic orbit.
		double speed = sqrt(2*mu/norm(x,y)) ;
		// The direction of velocity is calculated from the dy/dx but
		// it also need to be normalized. It is also set to face the 
		// origin.
		double velocity_direction_x = (x/abs(x))/ norm(1,x/2/R);
		double velocity_direction_y = (x/2/R   )/ norm(1,x/2/R);

		// Once the parameters are calculated, we can write them to
		// the ensemble data structure.
		s[b].mass() = mass_planet;
		s[b].x() = x  ,  s[b].vx() = -speed*velocity_direction_x ;
		s[b].y() = y  ,  s[b].vy() = -speed*velocity_direction_y ;
		s[b].z() = 0  ,  s[b].vz() = 0           ;
	}
	
	// End of the loop around systems.
}
// Now that the ensemble is created and filled, we save the results
// to the output file. Not that in real applications, we would want
// to integrate and examine the ensemble before writing it to a file
snapshot::save_text(ens,outputfn);
// On another note, we can also use snapshot::save function (with the same parameters)
// to write the ensemble to a binary file. Binary files are faster to read and
// write by C++ applications. The disadvantage is that they are not readable by
// Other applications. For manipulation of text and binary ensemble files
// see @ref SwarmExec.
}

// This concludes the tutorial. For complete listing of the file see @ref src/tutorials/tutorial_ensemble.cpp
