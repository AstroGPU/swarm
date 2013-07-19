/*************************************************************************
 * Copyright (C) 2011 by Saleh Dindar and the Swarm-NG Development Team  *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 3 of the License.        *
 *                                                                       *
 * This program is distributed in the hope t`hat it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program; if not, write to the                         *
 * Free Software Foundation, Inc.,                                       *
 * 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ************************************************************************/

/*! \file parabolic_collision.cpp
 *   \brief Create test case for collision
 * 
 * Put two planets on a collision course in a parabolic orbit
 * 
 * If we place the star at (0,0,0) then a parabola in XY plane will
 * be of the equation y = x^2/4r -r where r is the minimum
 * distance of the planet to the star on Y axis.
 *
 * the tangent line is (1, x/2r) which should be normalized.
 * 
 * The planets are put at x=-2 and x=+2. on the parabola and
 * they are set to meet at the bottom of the parabola. They meet
 * at (-1,0) at t=1.8856 .
 * 
 */

#include "swarm/swarm.h"
#include <iostream>
using namespace swarm;
using namespace std;

const int nbod = 3;
const int nsys = 64;

const double mass_planet = 0.0001;

/// For first two, the collision should be at 1.8856
/// at 3.77125 they end up at each others position.
const double xpos[4] = {  -2.0, +2.0, -4.0, +4.0 };


/// mu = GM , but we have set both G and M to 1.
const double mu = 1; 

/// the minimum distance of the planets to the star
/// on the parabolic orbit
const double R = 1; 

double norm(double x,double y){
	return sqrt(x*x+y*y);
}



int main(int argc, char* argv[]){
	if(argc <= 1){
		cout << "Usage: parabolic_collision <outputfilename>" << endl;
	}
	const string outputfn = argv[1];
	init(config());
	defaultEnsemble ens = defaultEnsemble::create(nbod,nsys);
	
	for(int i = 0; i < nsys ; i++){
		ensemble::SystemRef s = ens[i];
		s.id() = 0;
		s.time() = 0;
		s.set_active();
		
		// Stationary star at origin
		s[0].mass() = 1;
		s[0][0].pos() = 0, s[0][1].pos() = 0, s[0][2].pos() = 0 ;
		s[0][0].vel() = 0, s[0][1].vel() = 0, s[0][2].vel() = 0 ;
		
		
		// Planets on the parabola
		for(int b = 1; b < nbod; b++){
			s[b].mass() = mass_planet;
			
			double x = xpos[b-1], y = x*x/4/R - R ;
			double vmag = sqrt(2*mu/norm(x,y)) ;
			double vdirx = (x/abs(x))/norm(1,x/2/R), vdiry = abs(x)/2/R / norm(1,x/2/R);
			
			
			s[b][0].pos() = x  ,  s[b][0].vel() = -vmag*vdirx ;
			s[b][1].pos() = y  ,  s[b][1].vel() = -vmag*vdiry ;
			s[b][2].pos() = 0  ,  s[b][2].vel() = 0           ;
		}
		

	}
	snapshot::save_text(ens,outputfn);
}
