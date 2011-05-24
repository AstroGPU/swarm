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

#include "hp.hpp"
#include "static_accjerk.hpp"
#include "swarmlog.h"


namespace swarm {
namespace hp {

template<int nbod>
__global__ static void rkck_kernel(ensemble* ens,gpulog::device_log* dlog,double destination_time, double time_step){

	/////////////// FETCH LOCAL VARIABLES ///////////////////

	int thr = thread_in_system();

	if(sysid() >= ens->nsys()) { return; }
	ensemble::systemref sys ( (*ens)[sysid()] );

	// Body/Component Grid
	// Body number
	int b = thr / 3 ;
	// Component number
	int c = thr % 3 ;
	bool body_component_grid = b < nbod;

	// i,j pairs Grid
	int ij = thr;

	// shared memory allocation
	extern __shared__ char shared_mem[];

	double t_start = sys.time(), t = t_start;
	double t_end = min(t_start + destination_time,sys.time_end());

	// local information per component per body
	double pos = 0, vel = 0 , acc = 0, jerk = 0;
	if( body_component_grid )
		pos = sys[b].p(c), vel = sys[b].v(c);


	////////// INTEGRATION //////////////////////

	// Calculate acceleration and jerk
	Gravitation<nbod> calcForces(sys,shared_mem);
	calcForces(ij,b,c,pos,vel,acc,jerk);

	while(t < t_end){
		double h = min(time_step, t_end - t);

		const double b2 = 1.0 / 5.0;
		const double b3[] = { 3.0 / 40.0, 9.0 / 40.0 };
		const double b4[] = { 0.3, -0.9, 1.2 };
		const double b5[] = { -11.0 / 54.0, 2.5, -70.0 / 27.0, 35.0 / 27.0 };
		const double b6[] = { 1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0 };
		const double cff[]  = { 37.0 / 378.0, 0, 250.0 / 621.0, 125.0 / 594.0, 0 , 512.0 / 1771.0 } ;

		//// RKCK   integrate system
		double p0 = pos, v0 = vel;

		// Step 1
		double k1_acc = calcForces.acc(ij,b,c,p0,v0);
		double k1_vel = v0;

		double p1 = pos + h * b2 * k1_vel;
		double v1 = vel + h * b2 * k1_acc;

		// Step 2
		double k2_acc = calcForces.acc(ij,b,c,p1,v1);
		double k2_vel = v1;

		double p2 = pos + h * ( b3[0] * k1_vel + b3[1] * k2_vel );
		double v2 = vel + h * ( b3[0] * k1_acc + b3[1] * k2_acc );

		// Step 3
		double k3_acc = calcForces.acc(ij,b,c,p2,v2);
		double k3_vel = v2;

		double p3 = pos + h * ( b4[0] * k1_vel + b4[1] * k2_vel + b4[2] * k3_vel );
		double v3 = vel + h * ( b4[0] * k1_acc + b4[1] * k2_acc + b4[2] * k3_acc );

		// Step 4
		double k4_acc = calcForces.acc(ij,b,c,p3,v3);
		double k4_vel = v3;

		double p4 = pos + h * ( b5[0] * k1_vel + b5[1] * k2_vel + b5[2] * k3_vel + b5[3] * k4_vel );
		double v4 = vel + h * ( b5[0] * k1_acc + b5[1] * k2_acc + b5[2] * k3_acc + b5[3] * k4_acc );

		// Step 5
		double k5_acc = calcForces.acc(ij,b,c,p4,v4);
		double k5_vel = v4;

		double p5 = pos + h * ( b6[0] * k1_vel + b6[1] * k2_vel + b6[2] * k3_vel + b6[3] * k4_vel + b6[4] * k5_vel );
		double v5 = vel + h * ( b6[0] * k1_acc + b6[1] * k2_acc + b6[2] * k3_acc + b6[3] * k4_acc + b6[4] * k5_acc );

		// Step 6
		double k6_acc = calcForces.acc(ij,b,c,p5,v5);
		double k6_vel = v5;

		double p6 = pos + h * ( cff[0] * k1_vel + cff[1] * k2_vel + cff[2] * k3_vel + cff[3] * k4_vel + cff[4] * k5_vel + cff[5] * k6_vel );
		double v6 = vel + h * ( cff[0] * k1_acc + cff[1] * k2_acc + cff[2] * k3_acc + cff[3] * k4_acc + cff[4] * k5_acc + cff[5] * k6_acc );

		// Set the new positions and velocities
		pos = p6;
		vel = v6;

		t += h;

		if( body_component_grid )
			sys[b].p(c) = pos, sys[b].v(c) = vel;

		if(thr == 0) 
			if(log::needs_output(*ens, t, sysid()))
			{
				sys.set_time(t);
				log::output_system(*dlog, *ens, t, sysid());
			}

	}

	if(thr == 0) 
		sys.set_time(t);
}

class rkck: public integrator {
	private:
	double _time_step;

	public:
	rkck(const config& cfg): integrator(cfg),_time_step(0.001) {
		if(!cfg.count("time step")) ERROR("Integrator gpu_rkck requires a timestep ('time step' keyword in the config file).");
		_time_step = atof(cfg.at("time step").c_str());
	}

	template<int nbod>
	void launch_template(const double& destination_time)
	{
		if(_ens->nbod() == nbod) 
			rkck_kernel<nbod><<<gridDim(), threadDim(), shmemSize() >>>(_gpu_ens,_gpu_log,destination_time, _time_step);

	}

	virtual void launch_integrator(const double& destination_time){
			// flush CPU/GPU output logs
			log::flush(log::memory | log::if_full);

			if(_ens->nbod() <= 3){
				launch_template<3>(destination_time);
			} else {
				// How do we get an error message out of here?
				ERROR("Invalid number of bodies. (only up to 10 bodies per system)");
			}

			// flush CPU/GPU output logs
			log::flush(log::memory);
	}
};

/*!
 * \brief Factory to create double/single/mixed rkck gpu integrator based on precision
 *
 * @param[in] cfg configuration class
 *
 * @return        pointer to integrator cast to integrator*
 */
extern "C" integrator *create_hp_rkck(const config &cfg)
{
	return new rkck(cfg);
}

}
}

