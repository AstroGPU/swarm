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
#include "swarm/swarmplugin.h"
#include "keplerian.hpp"
#include "monitors/log_time_interval.hpp"

namespace swarm {

namespace gpu {
namespace bppt {

struct MVSPropagatorParams {
	double time_step;
	MVSPropagatorParams(const config& cfg){
		time_step = cfg.require("time_step", 0.0);
	}
};

//template<class T, class GravClass>
template<class T>
struct MVSPropagator {
	typedef MVSPropagatorParams params;
	static const int nbod = T::n;

	params _params;


	// Runtime variables
	ensemble::SystemRef& sys;
	Gravitation<T::n>& calcForces;
//	GravitationAccOnly<T::n>& calcForces;
//	GravClass& calcForces;	
	int b;
	int c;
	int ij;
	bool body_component_grid;
	bool first_thread_in_system;
	double sqrtGM;
	double max_timestep;


	GPUAPI MVSPropagator(const params& p,ensemble::SystemRef& s,
			Gravitation<T::n>& calc)
//			GravitationAccOnly<T::n>& calc)
//			GravClass calc)
		:_params(p),sys(s),calcForces(calc){}

	/// Shift into funky coordinate system (see A. Quillen's qymsym's tobary)
	GPUAPI void init()  { 
		sqrtGM = sqrt(sys[0].mass());

		if( body_component_grid )
		{
			double sump = 0., sumv = 0., mtot = 0.;

			// Find Center of mass and momentum
			for(int j=0;j<nbod;++j) {
				const double mj = sys[j].mass();
				mtot += mj;
				sump += mj*sys[j][c].pos();
				sumv += mj*sys[j][c].vel();
			}

			if(b==0) // For sun
				sys[b][c].vel() = sumv/mtot, 
					sys[b][c].pos() = sump/mtot;
			else     // For planets
				sys[b][c].vel() -= sumv/mtot,
					sys[b][c].pos() -= sys[0][c].pos();
		}

	}

	/// Shift back from funky coordinate system (see A. Quillen's qymsym's tobary)
	GPUAPI void shutdown() { 
		if ( body_component_grid ) {

			double sump = 0., sumv = 0., mtot;
			double m0, pc0, vc0;
			m0 = sys[0].mass();
			pc0 = sys[0][c].pos();
			vc0 = sys[0][c].vel();
			mtot = m0;
			for(int j=1;j<nbod;++j)
			{
				const double mj = sys[j].mass();
				mtot += mj;
				sump += mj*sys[j][c].pos();
				sumv += mj*sys[j][c].vel();
			}


			if(b==0) // For sun only
				sys[b][c].pos() -= sump/mtot,
					sys[b][c].vel() -= sumv/m0;
			else
				sys[b][c].pos() += pc0 - sump/mtot,
					sys[b][c].vel() += vc0;

		}

	}

	GPUAPI void drift_step(const double hby2) {
		double mv = 0;
		for(int j=1;j<nbod;++j)
			mv += sys[j].mass() * sys[j][c].vel();
		sys[b][c].pos() += mv*hby2/sys[0].mass();
	}


	GPUAPI void advance(){

		double H = min( max_timestep ,  _params.time_step );
		double hby2 = 0.5 * H;

		double acc = calcForces.acc_planets(ij,b,c);

		// Only operate on planets and not sun
		if( body_component_grid && b != 0 ) {


			// Step 1
			drift_step(hby2);

			// Step 2: Kick Step
			sys[b][c].vel() += hby2 * acc;

			// 3: Kepler Drift Step (Keplerian orbit about sun/central body)
			drift_kepler( sys[b][0].pos()
					     ,sys[b][1].pos()
						 ,sys[b][2].pos()

					     ,sys[b][0].vel()
					     ,sys[b][1].vel()
						 ,sys[b][2].vel()

						 ,sqrtGM, H
					);

			// TODO: check for close encounters here

			// Step 4: Kick Step
			sys[b][c].vel() += hby2 * acc;

			// Step 5
			drift_step(hby2);
		}

		if( first_thread_in_system ) 
			sys.time() += H;
	}
};

typedef gpulog::device_log L;
using namespace monitors;

integrator_plugin_initializer< generic< MVSPropagator, stop_on_ejection<L> > >
	mvs_prop_plugin("mvs"
			,"This is the integrator based on mvs propagator");


}
}
}

