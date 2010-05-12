/*************************************************************************
 * Copyright (C) 2010 by Eric Ford    and the Swarm-NG Development Team  *
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

/*! \file rk8pd.cu
 *  \brief declares prop_rk8pd for use with gpu_generic_integrator
 * 
 *  based on GSL's rk8pd.c and rk8pd_apply
*/

#include "swarm.h"
#include "rk8pd.h"
#include "swarmlog.h"

/// namespace for Swarm-NG library
namespace swarm {

// If want to declare __shared__ do it here (just like const variables for ensembles)

/*!  
 *  \brief propagator class for RKCK integrator on GPU: Advance the system by one time step.
 *
 *  CPU state and interface. Will be instantiated on construction of gpu_generic_integrator object. 
 *  Keep any data that need to reside on the CPU here.
 */
struct prop_rk8pd
{
	/*! 
	 * \brief GPU state and interface (per-grid). Will be passed as an argument to integration kernel. 
         *
	 * Any per-block read-only variables should be members of this structure.
         */
	struct gpu_t
	{
		//! per-block variables, time step
		double h;
		//! per-block variables, 
//		cuxDevicePtr<double, 3> k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;
		cuxDevicePtr<double, 4> k;
		//! per-block variables, initial & temporary coordinates, error term
//		cuxDevicePtr<double, 3> y0;
		cuxDevicePtr<double, 3> ytmp;
//		cuxDevicePtr<double, 3> yerr;

		// Prince-Dormand constants From GSL
		static const int order = 8;
		/*!
		 * \brief  GPU per-thread state and interface. Will be instantiated in the integration kernel. 
                 *
		 * Any per-thread variables should be members of this structure.
		 */
		struct thread_state_t
		{
			thread_state_t(const gpu_t &H, ensemble &ens, const int sys, double T, double Tend)
			{ }
		};

		/*! \brief compute accelerations from a temporary state in y
		 * @param[in]  ens   ensemble 
		 * @param[in]  y   temporary coordinates (3d array)
		 * @param[in]  sys    system id
 		 * @param[out] dydx  derivatives (3d array)
		 */
		__device__ void compute_acc(const ensemble &ens, const cuxDevicePtr<double, 3> y, const int sys, const int kidx)
		{
		typedef ThreeVector<double> V3;

		for ( unsigned int i=0;i<ens.nbod();++i )
		    {
		    V3 xi( y( sys,i,0 ), y( sys,i,1 ), y( sys,i,2 ) );
		    V3 ai(0.);
		    for (int j=ens.nbod()-1; j >= 0; j--)
			{
			if ( j==i ) continue; // Ignore body interacting with itself

			V3 dx(y(sys,j,0), y(sys,j,1), y(sys,j,2));  dx -= xi;
			double r2 = dx.MagnitudeSquared();
			double rinv = rsqrt ( r2 );
			rinv *= ens.mass ( sys,j );
			double rinv3 = rinv/r2;

			dx *= rinv3;
			ai += dx;
			}
		    k( kidx, sys, i, 0 ) = y( sys, i, 3 );
		    k( kidx, sys, i, 1 ) = y( sys, i, 4 );
		    k( kidx, sys, i, 2 ) = y( sys, i, 5 );
		    k( kidx, sys, i, 3 ) = ai.X();
		    k( kidx, sys, i, 4 ) = ai.Y();
		    k( kidx, sys, i, 5 ) = ai.Z();
		    } // end loop over bodies
		}

		/*! \brief compute accelerations from state in ensemble
		 * @param[in]  ens   ensemble 
		 * @param[in]  y     temporary coordinates (3d array)
		 * @param[in]  sys   system id
 		 * @param[out] dydx  derivatives (3d array)
		 */
		__device__ void compute_acc(const ensemble &ens, const int sys, const int kidx)
		{
		typedef ThreeVector<double> V3;

		for ( unsigned int i=0;i<ens.nbod();++i )
		    {
		    V3 xi( ens.x ( sys,i ), ens.y ( sys,i ), ens.z ( sys,i ) );
		    V3 ai(0.);
		    for (int j=ens.nbod()-1; j >= 0; j--)
			{
			if ( j==i ) continue; // Ignore body interacting with itself

			V3 dx(ens.x(sys,j), ens.y(sys,j), ens.z(sys,j));  dx -= xi;
			double r2 = dx.MagnitudeSquared();
			double rinv = rsqrt ( r2 );
			rinv *= ens.mass ( sys,j );
			double rinv3 = rinv/r2;

			dx *= rinv3;
			ai += dx;
			}
		    k( kidx, sys, i, 0 ) = ens.vx( sys, i);
		    k( kidx, sys, i, 1 ) = ens.vy( sys, i);
		    k( kidx, sys, i, 2 ) = ens.vz( sys, i);
		    k( kidx, sys, i, 3 ) = ai.X();
		    k( kidx, sys, i, 4 ) = ai.Y();
		    k( kidx, sys, i, 5 ) = ai.Z();
		    } // end loop over bodies
		}

		/*!
                 *  \brief Advance the system - this function must advance the system sys by one timestep, making sure that T does not exceed Tend.
		 *
		 * This function MUST return the new time of the system.
		 * This function MUST also call stop.test_body() for every body
		 * in the system, after that body has been advanced by a timestep.
                 *
  		 * see RKCK implementation GSL 
                 *
		 * @tparam stop_t ...
		 * @param[in,out] ens ensemble
		 * @param[in,out] pt ...
		 * @param[in] sys system ID
		 * @param[in] T start time
		 * @param[in] Tend destination time
		 * @param[in] stop ...
		 * @param[in] stop_ts ...
		 * @param[in] step  ...
		 * @return new time of the system
		 */
		template<typename stop_t>
		__device__ double advance(ensemble &ens, thread_state_t &pt, int sys, double T, double Tend, stop_t &stop, typename stop_t::thread_state_t &stop_ts, int step)
		{
//		const double ah[] = { 1.0 / 18.0, 1.0 / 12.0, 1.0 / 8.0, 5.0 / 16.0, 3.0 / 8.0, 59.0 / 400.0, 93.0 / 200.0, 5490023248.0 / 9719169821.0, 13.0 / 20.0, 1201146811.0 / 1299019798.0 };

			if(T >= Tend) { return T; }
			double hh = T + this->h <= Tend ? this->h : Tend - T;

			// k1 step
		        {
			compute_acc(ens,sys,0);       // at T
		        const double b21 = 1.0 / 18.0;
			double htmp = b21*hh;
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + htmp* k(0, sys, bod, 0);
			ytmp(sys, bod, 1) = ens.y (sys, bod) + htmp* k(0, sys, bod, 1);
			ytmp(sys, bod, 2) = ens.z (sys, bod) + htmp* k(0, sys, bod, 2);
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + htmp* k(0, sys, bod, 3);
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + htmp* k(0, sys, bod, 4);
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + htmp* k(0, sys, bod, 5);
			}
			} 
			// k2 step
			{
			const double b3[] = { 3.0 / 48.0, 1.0 / 16.0 };
			compute_acc(ens,ytmp,sys,1);  // at T + ah[0]*hh
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + hh* (b3[0]*k(0,sys, bod, 0) + b3[1]*k(1, sys, bod, 0));
			ytmp(sys, bod, 1) = ens.y (sys, bod) + hh* (b3[0]*k(0,sys, bod, 1) + b3[1]*k(1, sys, bod, 1));
			ytmp(sys, bod, 2) = ens.z (sys, bod) + hh* (b3[0]*k(0,sys, bod, 2) + b3[1]*k(1, sys, bod, 2));
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + hh* (b3[0]*k(0,sys, bod, 3) + b3[1]*k(1, sys, bod, 3));
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + hh* (b3[0]*k(0,sys, bod, 4) + b3[1]*k(1, sys, bod, 4));
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + hh* (b3[0]*k(0,sys, bod, 5) + b3[1]*k(1, sys, bod, 5));
			}
			}
			// k3 step
			{
			const double b4[] = { 1.0 / 32.0, 0.0, 3.0 / 32.0 };
			compute_acc(ens,ytmp,sys,2);  // at T + ah[1]*hh
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + hh* (b4[0]*k(0,sys, bod, 0) + b4[2]*k(2, sys, bod, 0) );
			ytmp(sys, bod, 1) = ens.y (sys, bod) + hh* (b4[0]*k(0,sys, bod, 1) + b4[2]*k(2, sys, bod, 1) );
			ytmp(sys, bod, 2) = ens.z (sys, bod) + hh* (b4[0]*k(0,sys, bod, 2) + b4[2]*k(2, sys, bod, 2) );
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + hh* (b4[0]*k(0,sys, bod, 3) + b4[2]*k(2, sys, bod, 3) );
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + hh* (b4[0]*k(0,sys, bod, 4) + b4[2]*k(2, sys, bod, 4) );
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + hh* (b4[0]*k(0,sys, bod, 5) + b4[2]*k(2, sys, bod, 5) );
			}
			}
			// k4 step
			{
			const double b5[] = { 5.0 / 16.0, 0.0, -75.0 / 64.0, 75.0 / 64.0 }; 
			compute_acc(ens,ytmp,sys,3);// at T + ah[2]*hh
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + hh* (b5[0]*k(0,sys, bod, 0) + b5[2]*k(2, sys, bod, 0) + b5[3]*k(3, sys, bod, 0) );
			ytmp(sys, bod, 1) = ens.y (sys, bod) + hh* (b5[0]*k(0,sys, bod, 1) + b5[2]*k(2, sys, bod, 1) + b5[3]*k(3, sys, bod, 1) );
			ytmp(sys, bod, 2) = ens.z (sys, bod) + hh* (b5[0]*k(0,sys, bod, 2) + b5[2]*k(2, sys, bod, 2) + b5[3]*k(3, sys, bod, 2) );
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + hh* (b5[0]*k(0,sys, bod, 3) + b5[2]*k(2, sys, bod, 3) + b5[3]*k(3, sys, bod, 3) );
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + hh* (b5[0]*k(0,sys, bod, 4) + b5[2]*k(2, sys, bod, 4) + b5[3]*k(3, sys, bod, 4) );
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + hh* (b5[0]*k(0,sys, bod, 5) + b5[2]*k(2, sys, bod, 5) + b5[3]*k(3, sys, bod, 5) );
			}
			}
			// k5 step
			{
			const double b6[] = { 3.0 / 80.0, 0.0, 0.0, 3.0 / 16.0, 3.0 / 20.0 };
			compute_acc(ens,ytmp,sys,4);// at T + ah[3]*hh
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + hh* (b6[0]*k(0,sys, bod, 0) + b6[3]*k(3, sys, bod, 0) + b6[4]*k(4, sys, bod, 0));
			ytmp(sys, bod, 1) = ens.y (sys, bod) + hh* (b6[0]*k(0,sys, bod, 1) + b6[3]*k(3, sys, bod, 1) + b6[4]*k(4, sys, bod, 1));
			ytmp(sys, bod, 2) = ens.z (sys, bod) + hh* (b6[0]*k(0,sys, bod, 2) + b6[3]*k(3, sys, bod, 2) + b6[4]*k(4, sys, bod, 2));
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + hh* (b6[0]*k(0,sys, bod, 3) + b6[3]*k(3, sys, bod, 3) + b6[4]*k(4, sys, bod, 3));
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + hh* (b6[0]*k(0,sys, bod, 4) + b6[3]*k(3, sys, bod, 4) + b6[4]*k(4, sys, bod, 4));
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + hh* (b6[0]*k(0,sys, bod, 5) + b6[3]*k(3, sys, bod, 5) + b6[4]*k(4, sys, bod, 5));
			}
			}
			// k6 step
			{
			const double b7[] = { 29443841.0 / 614563906.0, 0.0, 0.0, 77736538.0 / 692538347.0, -28693883.0 / 1125000000.0, 23124283.0 / 1800000000.0};
 			compute_acc(ens,ytmp,sys,5);// at T + ah[4]*hh
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + hh* (b7[0]*k(0,sys, bod, 0) + b7[3]*k(3, sys, bod, 0) + b7[4]*k(4, sys, bod, 0) + b7[5]*k(, sys, bod, 0));
			ytmp(sys, bod, 1) = ens.y (sys, bod) + hh* (b7[0]*k(0,sys, bod, 1) + b7[3]*k(3, sys, bod, 1) + b7[4]*k(4, sys, bod, 1) + b7[5]*k(, sys, bod, 1));
			ytmp(sys, bod, 2) = ens.z (sys, bod) + hh* (b7[0]*k(0,sys, bod, 2) + b7[3]*k(3, sys, bod, 2) + b7[4]*k(4, sys, bod, 2) + b7[5]*k(, sys, bod, 2));
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + hh* (b7[0]*k(0,sys, bod, 3) + b7[3]*k(3, sys, bod, 3) + b7[4]*k(4, sys, bod, 3) + b7[5]*k(, sys, bod, 3));
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + hh* (b7[0]*k(0,sys, bod, 4) + b7[3]*k(3, sys, bod, 4) + b7[4]*k(4, sys, bod, 4) + b7[5]*k(, sys, bod, 4));
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + hh* (b7[0]*k(0,sys, bod, 5) + b7[3]*k(3, sys, bod, 5) + b7[4]*k(4, sys, bod, 5) + b7[5]*k(, sys, bod, 5));
			}
			}
			// k7 step
			{
			const double b8[] = { 16016141.0 / 946692911.0, 0.0, 0.0, 61564180.0 / 158732637.0, 22789713.0 / 633445777.0, 545815736.0 / 2771057229.0, -180193667.0 / 1043307555.0}; 
			compute_acc(ens,ytmp,sys,6);// at T + ah[5]*hh
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + hh* (b8[0]*k(0,sys, bod, 0) + b8[3]*k(3, sys, bod, 0) + b8[4]*k(4, sys, bod, 0) + b8[5]*k(, sys, bod, 0) + b8[6]*k(6, sys, bod, 0));
			ytmp(sys, bod, 1) = ens.y (sys, bod) + hh* (b8[0]*k(0,sys, bod, 1) + b8[3]*k(3, sys, bod, 1) + b8[4]*k(4, sys, bod, 1) + b8[5]*k(, sys, bod, 1) + b8[6]*k(6, sys, bod, 1));
			ytmp(sys, bod, 2) = ens.z (sys, bod) + hh* (b8[0]*k(0,sys, bod, 2) + b8[3]*k(3, sys, bod, 2) + b8[4]*k(4, sys, bod, 2) + b8[5]*k(, sys, bod, 2) + b8[6]*k(6, sys, bod, 2));
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + hh* (b8[0]*k(0,sys, bod, 3) + b8[3]*k(3, sys, bod, 3) + b8[4]*k(4, sys, bod, 3) + b8[5]*k(, sys, bod, 3) + b8[6]*k(6, sys, bod, 3));
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + hh* (b8[0]*k(0,sys, bod, 4) + b8[3]*k(3, sys, bod, 4) + b8[4]*k(4, sys, bod, 4) + b8[5]*k(, sys, bod, 4) + b8[6]*k(6, sys, bod, 4));
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + hh* (b8[0]*k(0,sys, bod, 5) + b8[3]*k(3, sys, bod, 5) + b8[4]*k(4, sys, bod, 5) + b8[5]*k(, sys, bod, 5) + b8[6]*k(6, sys, bod, 5));
			}
			}
			// k8 step
			{
			const double b9[] = { 39632708.0 / 573591083.0, 0.0, 0.0, -433636366.0 / 683701615.0, -421739975.0 / 2616292301.0, 100302831.0 / 723423059.0, 790204164.0 / 839813087.0, 800635310.0 / 3783071287.0};
			compute_acc(ens,ytmp,sys,7);// at T + ah[6]*hh
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + hh* (b9[0]*k(0,sys, bod, 0) + b9[3]*k(3, sys, bod, 0) + b9[4]*k(4, sys, bod, 0) + b9[5]*k(, sys, bod, 0) + b9[6]*k(6, sys, bod, 0) + b9[7]*k(7, sys, bod, 0));
			ytmp(sys, bod, 1) = ens.y (sys, bod) + hh* (b9[0]*k(0,sys, bod, 1) + b9[3]*k(3, sys, bod, 1) + b9[4]*k(4, sys, bod, 1) + b9[5]*k(, sys, bod, 1) + b9[6]*k(6, sys, bod, 1) + b9[7]*k(7, sys, bod, 1));
			ytmp(sys, bod, 2) = ens.z (sys, bod) + hh* (b9[0]*k(0,sys, bod, 2) + b9[3]*k(3, sys, bod, 2) + b9[4]*k(4, sys, bod, 2) + b9[5]*k(, sys, bod, 2) + b9[6]*k(6, sys, bod, 2) + b9[7]*k(7, sys, bod, 2));
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + hh* (b9[0]*k(0,sys, bod, 3) + b9[3]*k(3, sys, bod, 3) + b9[4]*k(4, sys, bod, 3) + b9[5]*k(, sys, bod, 3) + b9[6]*k(6, sys, bod, 3) + b9[7]*k(7, sys, bod, 3));
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + hh* (b9[0]*k(0,sys, bod, 4) + b9[3]*k(3, sys, bod, 4) + b9[4]*k(4, sys, bod, 4) + b9[5]*k(, sys, bod, 4) + b9[6]*k(6, sys, bod, 4) + b9[7]*k(7, sys, bod, 4));
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + hh* (b9[0]*k(0,sys, bod, 5) + b9[3]*k(3, sys, bod, 5) + b9[4]*k(4, sys, bod, 5) + b9[5]*k(, sys, bod, 5) + b9[6]*k(6, sys, bod, 5) + b9[7]*k(7, sys, bod, 5));
			}
			}
			// k9 step
			{
			const double b10[] = {  246121993.0 / 1340847787.0, 0.0, 0.0, -37695042795.0 / 15268766246.0, -309121744.0 / 1061227803.0, -12992083.0 / 490766935.0, 6005943493.0 / 2108947869.0, 393006217.0 / 1396673457.0, 123872331.0 / 1001029789.0 };
			compute_acc(ens,ytmp,sys,8);// at T + ah[7]*hh
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + hh* (b10[0]*k(0,sys, bod, 0) + b10[3]*k(3, sys, bod, 0) + b10[4]*k(4, sys, bod, 0) + b10[5]*k(, sys, bod, 0) + b10[6]*k(6, sys, bod, 0) + b10[7]*k(7, sys, bod, 0) + b10[8]*k(8, sys, bod, 0));
			ytmp(sys, bod, 1) = ens.y (sys, bod) + hh* (b10[0]*k(0,sys, bod, 1) + b10[3]*k(3, sys, bod, 1) + b10[4]*k(4, sys, bod, 1) + b10[5]*k(, sys, bod, 1) + b10[6]*k(6, sys, bod, 1) + b10[7]*k(7, sys, bod, 1) + b10[8]*k(8, sys, bod, 1));
			ytmp(sys, bod, 2) = ens.z (sys, bod) + hh* (b10[0]*k(0,sys, bod, 2) + b10[3]*k(3, sys, bod, 2) + b10[4]*k(4, sys, bod, 2) + b10[5]*k(, sys, bod, 2) + b10[6]*k(6, sys, bod, 2) + b10[7]*k(7, sys, bod, 2) + b10[8]*k(8, sys, bod, 2));
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + hh* (b10[0]*k(0,sys, bod, 3) + b10[3]*k(3, sys, bod, 3) + b10[4]*k(4, sys, bod, 3) + b10[5]*k(, sys, bod, 3) + b10[6]*k(6, sys, bod, 3) + b10[7]*k(7, sys, bod, 3) + b10[8]*k(8, sys, bod, 3));
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + hh* (b10[0]*k(0,sys, bod, 4) + b10[3]*k(3, sys, bod, 4) + b10[4]*k(4, sys, bod, 4) + b10[5]*k(, sys, bod, 4) + b10[6]*k(6, sys, bod, 4) + b10[7]*k(7, sys, bod, 4) + b10[8]*k(8, sys, bod, 4));
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + hh* (b10[0]*k(0,sys, bod, 5) + b10[3]*k(3, sys, bod, 5) + b10[4]*k(4, sys, bod, 5) + b10[5]*k(, sys, bod, 5) + b10[6]*k(6, sys, bod, 5) + b10[7]*k(7, sys, bod, 5) + b10[8]*k(8, sys, bod, 5));
			}
			}
			// k10 step
			{
			const double b11[] = { -1028468189.0 / 846180014.0, 0.0, 0.0, 8478235783.0 / 508512852.0, 1311729495.0 / 1432422823.0, -10304129995.0 / 1701304382.0, -48777925059.0 / 3047939560.0, 15336726248.0 / 1032824649.0, -45442868181.0 / 3398467696.0, 3065993473.0 / 597172653.0 }; 
			compute_acc(ens,ytmp,sys,9);// at T + ah[8]*hh
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + hh* (b11[0]*k(0,sys, bod, 0) + b11[3]*k(3, sys, bod, 0) + b11[4]*k(4, sys, bod, 0) + b11[5]*k(, sys, bod, 0) + b11[6]*k(6, sys, bod, 0) + b11[7]*k(7, sys, bod, 0) + b11[8]*k(8, sys, bod, 0) + b11[9]*k(9, sys, bod, 0));
			ytmp(sys, bod, 1) = ens.y (sys, bod) + hh* (b11[0]*k(0,sys, bod, 1) + b11[3]*k(3, sys, bod, 1) + b11[4]*k(4, sys, bod, 1) + b11[5]*k(, sys, bod, 1) + b11[6]*k(6, sys, bod, 1) + b11[7]*k(7, sys, bod, 1) + b11[8]*k(8, sys, bod, 1) + b11[9]*k(9, sys, bod, 1));
			ytmp(sys, bod, 2) = ens.z (sys, bod) + hh* (b11[0]*k(0,sys, bod, 2) + b11[3]*k(3, sys, bod, 2) + b11[4]*k(4, sys, bod, 2) + b11[5]*k(, sys, bod, 2) + b11[6]*k(6, sys, bod, 2) + b11[7]*k(7, sys, bod, 2) + b11[8]*k(8, sys, bod, 2) + b11[9]*k(9, sys, bod, 2));
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + hh* (b11[0]*k(0,sys, bod, 3) + b11[3]*k(3, sys, bod, 3) + b11[4]*k(4, sys, bod, 3) + b11[5]*k(, sys, bod, 3) + b11[6]*k(6, sys, bod, 3) + b11[7]*k(7, sys, bod, 3) + b11[8]*k(8, sys, bod, 3) + b11[9]*k(9, sys, bod, 3));
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + hh* (b11[0]*k(0,sys, bod, 4) + b11[3]*k(3, sys, bod, 4) + b11[4]*k(4, sys, bod, 4) + b11[5]*k(, sys, bod, 4) + b11[6]*k(6, sys, bod, 4) + b11[7]*k(7, sys, bod, 4) + b11[8]*k(8, sys, bod, 4) + b11[9]*k(9, sys, bod, 4));
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + hh* (b11[0]*k(0,sys, bod, 5) + b11[3]*k(3, sys, bod, 5) + b11[4]*k(4, sys, bod, 5) + b11[5]*k(, sys, bod, 5) + b11[6]*k(6, sys, bod, 5) + b11[7]*k(7, sys, bod, 5) + b11[8]*k(8, sys, bod, 5) + b11[9]*k(9, sys, bod, 5));
			}
			}
			// k11 step
			{
			const double b12[] = { 185892177.0 / 718116043.0, 0.0, 0.0, -3185094517.0 / 667107341.0, -477755414.0 / 1098053517.0, -703635378.0 / 230739211.0, 5731566787.0 / 1027545527.0, 5232866602.0 / 850066563.0, -4093664535.0 / 808688257.0, 3962137247.0 / 1805957418.0, 65686358.0 / 487910083.0};
			compute_acc(ens,ytmp,sys,10);// at T + ah[9]*hh
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + hh* (b12[0]*k(0,sys, bod, 0) + b12[3]*k(3, sys, bod, 0) + b12[4]*k(4, sys, bod, 0) + b12[5]*k(, sys, bod, 0) + b12[6]*k(6, sys, bod, 0) + b12[7]*k(7, sys, bod, 0) + b12[8]*k(8, sys, bod, 0) + b12[9]*k(9, sys, bod, 0) + b12[10]*k(10, sys, bod, 0));
			ytmp(sys, bod, 1) = ens.y (sys, bod) + hh* (b12[0]*k(0,sys, bod, 1) + b12[3]*k(3, sys, bod, 1) + b12[4]*k(4, sys, bod, 1) + b12[5]*k(, sys, bod, 1) + b12[6]*k(6, sys, bod, 1) + b12[7]*k(7, sys, bod, 1) + b12[8]*k(8, sys, bod, 1) + b12[9]*k(9, sys, bod, 1) + b12[10]*k(10, sys, bod, 1));
			ytmp(sys, bod, 2) = ens.z (sys, bod) + hh* (b12[0]*k(0,sys, bod, 2) + b12[3]*k(3, sys, bod, 2) + b12[4]*k(4, sys, bod, 2) + b12[5]*k(, sys, bod, 2) + b12[6]*k(6, sys, bod, 2) + b12[7]*k(7, sys, bod, 2) + b12[8]*k(8, sys, bod, 2) + b12[9]*k(9, sys, bod, 2) + b12[10]*k(10, sys, bod, 2));
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + hh* (b12[0]*k(0,sys, bod, 3) + b12[3]*k(3, sys, bod, 3) + b12[4]*k(4, sys, bod, 3) + b12[5]*k(, sys, bod, 3) + b12[6]*k(6, sys, bod, 3) + b12[7]*k(7, sys, bod, 3) + b12[8]*k(8, sys, bod, 3) + b12[9]*k(9, sys, bod, 3) + b12[10]*k(10, sys, bod, 3));
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + hh* (b12[0]*k(0,sys, bod, 4) + b12[3]*k(3, sys, bod, 4) + b12[4]*k(4, sys, bod, 4) + b12[5]*k(, sys, bod, 4) + b12[6]*k(6, sys, bod, 4) + b12[7]*k(7, sys, bod, 4) + b12[8]*k(8, sys, bod, 4) + b12[9]*k(9, sys, bod, 4) + b12[10]*k(10, sys, bod, 4));
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + hh* (b12[0]*k(0,sys, bod, 5) + b12[3]*k(3, sys, bod, 5) + b12[4]*k(4, sys, bod, 5) + b12[5]*k(, sys, bod, 5) + b12[6]*k(6, sys, bod, 5) + b12[7]*k(7, sys, bod, 5) + b12[8]*k(8, sys, bod, 5) + b12[9]*k(9, sys, bod, 5) + b12[10]*k(10, sys, bod, 5));
			}
			}
			// k12 step
			{
			const double b13[] = {  403863854.0 / 491063109.0, 0.0, 0.0, -5068492393.0 / 434740067.0, -411421997.0 / 543043805.0, 652783627.0 / 914296604.0, 11173962825.0 / 925320556.0, -13158990841.0 / 6184727034.0, 3936647629.0 / 1978049680.0, -160528059.0 / 685178525.0, 248638103.0 / 1413531060.0, 0.0};
			compute_acc(ens,ytmp,sys,11);// at T + hh
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + hh* (b13[0]*k(0,sys, bod, 0) + b13[3]*k(3, sys, bod, 0) + b13[4]*k(4, sys, bod, 0) + b13[5]*k(, sys, bod, 0) + b13[6]*k(6, sys, bod, 0) + b13[7]*k(7, sys, bod, 0) + b13[8]*k(8, sys, bod, 0) + b13[9]*k(9, sys, bod, 0) + b13[10]*k(10, sys, bod, 0) + b13[11]*k(11, sys, bod, 0));
			ytmp(sys, bod, 1) = ens.y (sys, bod) + hh* (b13[0]*k(0,sys, bod, 1) + b13[3]*k(3, sys, bod, 1) + b13[4]*k(4, sys, bod, 1) + b13[5]*k(, sys, bod, 1) + b13[6]*k(6, sys, bod, 1) + b13[7]*k(7, sys, bod, 1) + b13[8]*k(8, sys, bod, 1) + b13[9]*k(9, sys, bod, 1) + b13[10]*k(10, sys, bod, 1) + b13[11]*k(11, sys, bod, 1));
			ytmp(sys, bod, 2) = ens.z (sys, bod) + hh* (b13[0]*k(0,sys, bod, 2) + b13[3]*k(3, sys, bod, 2) + b13[4]*k(4, sys, bod, 2) + b13[5]*k(, sys, bod, 2) + b13[6]*k(6, sys, bod, 2) + b13[7]*k(7, sys, bod, 2) + b13[8]*k(8, sys, bod, 2) + b13[9]*k(9, sys, bod, 2) + b13[10]*k(10, sys, bod, 2) + b13[11]*k(11, sys, bod, 2));
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + hh* (b13[0]*k(0,sys, bod, 3) + b13[3]*k(3, sys, bod, 3) + b13[4]*k(4, sys, bod, 3) + b13[5]*k(, sys, bod, 3) + b13[6]*k(6, sys, bod, 3) + b13[7]*k(7, sys, bod, 3) + b13[8]*k(8, sys, bod, 3) + b13[9]*k(9, sys, bod, 3) + b13[10]*k(10, sys, bod, 3) + b13[11]*k(11, sys, bod, 3));
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + hh* (b13[0]*k(0,sys, bod, 4) + b13[3]*k(3, sys, bod, 4) + b13[4]*k(4, sys, bod, 4) + b13[5]*k(, sys, bod, 4) + b13[6]*k(6, sys, bod, 4) + b13[7]*k(7, sys, bod, 4) + b13[8]*k(8, sys, bod, 4) + b13[9]*k(9, sys, bod, 4) + b13[10]*k(10, sys, bod, 4) + b13[11]*k(11, sys, bod, 4));
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + hh* (b13[0]*k(0,sys, bod, 5) + b13[3]*k(3, sys, bod, 5) + b13[4]*k(4, sys, bod, 5) + b13[5]*k(, sys, bod, 5) + b13[6]*k(6, sys, bod, 5) + b13[7]*k(7, sys, bod, 5) + b13[8]*k(8, sys, bod, 5) + b13[9]*k(9, sys, bod, 5) + b13[10]*k(10, sys, bod, 5) + b13[11]*k(11, sys, bod, 5));
			}
			}
			// k13 step & final sum
			{
			const double Abar[] = { 14005451.0 / 335480064.0, 0.0, 0.0, 0.0, 0.0, -59238493.0 / 1068277825.0, 181606767.0 / 758867731.0, 561292985.0 / 797845732.0, -1041891430.0 / 1371343529.0, 760417239.0 / 1151165299.0, 118820643.0 / 751138087.0, -528747749.0 / 2220607170.0, 1.0 / 4.0 };
//			const double A[] = { 13451932.0 / 455176623.0, 0.0, 0.0, 0.0, 0.0, -808719846.0 / 976000145.0, 1757004468.0 / 5645159321.0, 656045339.0 / 265891186.0, -3867574721.0 / 1518517206.0, 465885868.0 / 322736535.0, 53011238.0 / 667516719.0, 2.0 / 45.0 };
			compute_acc(ens,ytmp,sys,12);// at T + hh
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			   double ksum8;
//			   double ksum7;
			   ksum8 = Abar[0] * k(0,sys, bod, 0) + Abar[5] * k(, sys, bod, 0) + Abar[6] * k(6, sys, bod, 0) + Abar[7] * k(7, sys, bod, 0) + Abar[8] * k(8, sys, bod, 0) + Abar[9] * k(9, sys, bod, 0) + Abar[10] * k(10, sys, bod, 0) + Abar[11] * k(11, sys, bod, 0) + Abar[12] * k(12, sys, bod, 0);
//			   ksum7 = A[0] * k(0,sys, bod, 0) + A[5] * k(, sys, bod, 0) + A[6] * k(6, sys, bod, 0) + A[7] * k(7, sys, bod, 0) + A[8] * k(8, sys, bod, 0) + A[9] * k(9, sys, bod, 0) + A[10] * k(10, sys, bod, 0) + A[11] * k(11, sys, bod, 0);
			   ens.x(sys,bod) += hh * ksum8; 
			   // yerr(sys,bod, 0) = hh * (ksum7-ksum8);

			   ksum8 = Abar[0] * k(0,sys, bod, 1) + Abar[5] * k(, sys, bod, 1) + Abar[6] * k(6, sys, bod, 1) + Abar[7] * k(7, sys, bod, 1) + Abar[8] * k(8, sys, bod, 1) + Abar[9] * k(9, sys, bod, 1) + Abar[10] * k(10, sys, bod, 1) + Abar[11] * k(11, sys, bod, 1) + Abar[12] * k(12, sys, bod, 1);
//			   ksum7 = A[0] * k(0,sys, bod, 1) + A[5] * k(, sys, bod, 1) + A[6] * k(6, sys, bod, 1) + A[7] * k(7, sys, bod, 1) + A[8] * k(8, sys, bod, 1) + A[9] * k(9, sys, bod, 1) + A[10] * k(10, sys, bod, 1) + A[11] * k(11, sys, bod, 1);
			   ens.y(sys,bod) += hh * ksum8; 
			   // yerr(sys,bod, 1) = hh * (ksum7-ksum8);

			   ksum8 = Abar[0] * k(0,sys, bod, 2) + Abar[5] * k(, sys, bod, 2) + Abar[6] * k(6, sys, bod, 2) + Abar[7] * k(7, sys, bod, 2) + Abar[8] * k(8, sys, bod, 2) + Abar[9] * k(9, sys, bod, 2) + Abar[10] * k(10, sys, bod, 2) + Abar[11] * k(11, sys, bod, 2) + Abar[12] * k(12, sys, bod, 2);
//			   ksum7 = A[0] * k(0,sys, bod, 2) + A[5] * k(, sys, bod, 2) + A[6] * k(6, sys, bod, 2) + A[7] * k(7, sys, bod, 2) + A[8] * k(8, sys, bod, 2) + A[9] * k(9, sys, bod, 2) + A[10] * k(10, sys, bod, 2) + A[11] * k(11, sys, bod, 2);
			   ens.z(sys,bod) += hh * ksum8; 
			   // yerr(sys,bod, 2) = hh * (ksum7-ksum8);

			   ksum8 = Abar[0] * k(0,sys, bod, 3) + Abar[5] * k(, sys, bod, 3) + Abar[6] * k(6, sys, bod, 3) + Abar[7] * k(7, sys, bod, 3) + Abar[8] * k(8, sys, bod, 3) + Abar[9] * k(9, sys, bod, 3) + Abar[10] * k(10, sys, bod, 3) + Abar[11] * k(11, sys, bod, 3) + Abar[12] * k(12, sys, bod, 3);
//			   ksum7 = A[0] * k(0,sys, bod, 3) + A[5] * k(, sys, bod, 3) + A[6] * k(6, sys, bod, 3) + A[7] * k(7, sys, bod, 3) + A[8] * k(8, sys, bod, 3) + A[9] * k(9, sys, bod, 3) + A[10] * k(10, sys, bod, 3) + A[11] * k(11, sys, bod, 3);
			   ens.vx(sys,bod) += hh * ksum8; 
			   // yerr(sys,bod, 3) = hh * (ksum7-ksum8);

			   ksum8 = Abar[0] * k(0,sys, bod, 4) + Abar[5] * k(, sys, bod, 4) + Abar[6] * k(6, sys, bod, 4) + Abar[7] * k(7, sys, bod, 4) + Abar[8] * k(8, sys, bod, 4) + Abar[9] * k(9, sys, bod, 4) + Abar[10] * k(10, sys, bod, 4) + Abar[11] * k(11, sys, bod, 4) + Abar[12] * k(12, sys, bod, 4);
//			   ksum7 = A[0] * k(0,sys, bod, 4) + A[5] * k(, sys, bod, 4) + A[6] * k(6, sys, bod, 4) + A[7] * k(7, sys, bod, 4) + A[8] * k(8, sys, bod, 4) + A[9] * k(9, sys, bod, 4) + A[10] * k(10, sys, bod, 4) + A[11] * k(11, sys, bod, 4);
			   ens.vy(sys,bod) += hh * ksum8; 
			   // yerr(sys,bod, 4) = hh * (ksum7-ksum8);

			   ksum8 = Abar[0] * k(0,sys, bod, 5) + Abar[5] * k(, sys, bod, 5) + Abar[6] * k(6, sys, bod, 5) + Abar[7] * k(7, sys, bod, 5) + Abar[8] * k(8, sys, bod, 5) + Abar[9] * k(9, sys, bod, 5) + Abar[10] * k(10, sys, bod, 5) + Abar[11] * k(11, sys, bod, 5) + Abar[12] * k(12, sys, bod, 5);
//			   ksum7 = A[0] * k(0,sys, bod, 5) + A[5] * k(, sys, bod, 5) + A[6] * k(6, sys, bod, 5) + A[7] * k(7, sys, bod, 5) + A[8] * k(8, sys, bod, 5) + A[9] * k(9, sys, bod, 5) + A[10] * k(10, sys, bod, 5) + A[11] * k(11, sys, bod, 5);
			   ens.vz(sys,bod) += hh * ksum8; 
			   // yerr(sys,bod, 5) = hh * (ksum7-ksum8);

			stop.test_body(stop_ts, ens, sys, bod, T+hh, ens.x(sys,bod), ens.y(sys,bod), ens.z(sys,bod), ens.vx(sys,bod), ens.vy(sys,bod), ens.vz(sys,bod));			
			}
			}
			return T + hh;
		}

	};

	//! CPU state and interface
//	cuxDeviceAutoPtr<double, 3> k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;
	cuxDeviceAutoPtr<double, 4> k;
//	cuxDeviceAutoPtr<double, 3> y0;
	cuxDeviceAutoPtr<double, 3> ytmp;
//	cuxDeviceAutoPtr<double, 3> yerr;
	gpu_t gpu_obj;

	/*!
         * \brief initialize temporary variables for ensemble ens. 
         *
         * This function should initialize any temporary state that is needed for integration of ens. 
	 * It will be called from gpu_generic_integrator, but only if ens.last_integrator() != this. 
         * If any temporary state exists from previous invocation of this function, it should be deallocated and the new state (re)allocated.
	 */
	void initialize(ensemble &ens)
	{
		// Here you'd initialize the object to be passed to the kernel, or
		// upload any temporary data you need to constant/texture/global memory
		k.realloc(13, ens.nsys(), ens.nbod(), 6);
/*
		k1.realloc(ens.nsys(), ens.nbod(), 6);
		k2.realloc(ens.nsys(), ens.nbod(), 6);
		k3.realloc(ens.nsys(), ens.nbod(), 6);
		k4.realloc(ens.nsys(), ens.nbod(), 6);
		k5.realloc(ens.nsys(), ens.nbod(), 6);
		k6.realloc(ens.nsys(), ens.nbod(), 6);
		k7.realloc(ens.nsys(), ens.nbod(), 6);
		k8.realloc(ens.nsys(), ens.nbod(), 6);
		k9.realloc(ens.nsys(), ens.nbod(), 6);
		k10.realloc(ens.nsys(), ens.nbod(), 6);
		k11.realloc(ens.nsys(), ens.nbod(), 6);
		k12.realloc(ens.nsys(), ens.nbod(), 6);
		k13.realloc(ens.nsys(), ens.nbod(), 6);
*/
//		y0.realloc(ens.nsys(), ens.nbod(), 6);
		ytmp.realloc(ens.nsys(), ens.nbod(), 6);
//		yerr.realloc(ens.nsys(), ens.nbod(), 6);
		gpu_obj.k= k;
/*
		gpu_obj.k1= k1;
		gpu_obj.k2= k2;
		gpu_obj.k3= k3;
		gpu_obj.k4= k4;
		gpu_obj.k5= k5;
		gpu_obj.k6= k6;
		gpu_obj.k7= k7;
		gpu_obj.k8= k8;
		gpu_obj.k9= k9;
		gpu_obj.k10= k10;
		gpu_obj.k11= k11;
		gpu_obj.k12= k12;
		gpu_obj.k13= k13;
*/
//		gpu_obj.y0= y0;
		gpu_obj.ytmp= ytmp;
//		gpu_obj.yerr= yerr;
	}

	/*!
	 * \brief constructor 
         * 
         * Constructor will be passed the cfg object with the contents of integrator configuration file. 
         * It will be called during construction of gpu_generic_integrator. 
         * It should load any values necessary for initialization.
	 */
	prop_rk8pd(const config &cfg)
	{
		if(!cfg.count("time step")) ERROR("Integrator gpu_rk8pd requires a timestep ('time step' keyword in the config file).");
		gpu_obj.h = atof(cfg.at("time step").c_str());
	}

	/*!
         * \brief Cast operator for gpu_t.
         *
         * This operator must return the gpu_t object to be passed to integration kernel. 
         * It is called once per kernel invocation.
	 * @return gpu_t object to be passed to integration kernel.
	 */
	operator gpu_t()
	{
		return gpu_obj;
	}
};


/*!
 * \brief factory function to create an integrator 
 * 	  
 * This factory uses the gpu_generic_integrator class
 * with the propagator rk8pd and the stopper stop_on_ejection
 * 
 * @param[in] cfg contains configuration data for gpu_generic_integrator
 */
extern "C" integrator *create_gpu_rk8pd(const config &cfg)
{
	return new gpu_generic_integrator<stop_on_ejection, prop_rk8pd>(cfg);
       // return new gpu_generic_integrator< stop_on_crossing_orbit_or_close_approach, prop_rk8pd>(cfg);
}

} // end namespace swarm
