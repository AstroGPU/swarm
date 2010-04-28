/*! \file rkck.cu
 *  \brief declares prop_rkck for use with gpu_generic_integrator
 * 
 *  based on GSL's rkck.c and rkck_apply
*/

#include "swarm.h"
#include "rkck.h"
#include "swarmlog.h"

/// namespace for Swarm-NG library
namespace swarm {

/*!  
 *  \brief propagator class for RKCK integrator on GPU: Advance the system by one time step.
 *
 *  CPU state and interface. Will be instantiated on construction of gpu_generic_integrator object. 
 *  Keep any data that need to reside on the CPU here.
 */
struct prop_rkck
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
		cuxDevicePtr<double, 3> k1, k2, k3, k4, k5, k6;
		//! per-block variables, initial & temporary coordinates, error term
//		cuxDevicePtr<double, 3> y0;
		cuxDevicePtr<double, 3> ytmp;
//		cuxDevicePtr<double, 3> yerr;

		// Cash-Karp constants From GSL
		static const int order = 5;
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
		__device__ void compute_acc(ensemble &ens, cuxDevicePtr<double, 3> y, int sys, cuxDevicePtr<double, 3> dydx)
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
		    dydx ( sys, i, 0 ) = y( sys, i, 3 );
		    dydx ( sys, i, 1 ) = y( sys, i, 4 );
		    dydx ( sys, i, 2 ) = y( sys, i, 5 );
		    dydx ( sys, i, 3 ) = ai.X();
		    dydx ( sys, i, 4 ) = ai.Y();
		    dydx ( sys, i, 5 ) = ai.Z();
		    } // end loop over bodies
		}

		/*! \brief compute accelerations from state in ensemble
		 * @param[in]  ens   ensemble 
		 * @param[in]  y     temporary coordinates (3d array)
		 * @param[in]  sys   system id
 		 * @param[out] dydx  derivatives (3d array)
		 */
		__device__ void compute_acc(ensemble &ens, int sys, cuxDevicePtr<double, 3> dydx)
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
		    dydx ( sys, i, 0 ) = ens.vx( sys, i);
		    dydx ( sys, i, 1 ) = ens.vy( sys, i);
		    dydx ( sys, i, 2 ) = ens.vz( sys, i);
		    dydx ( sys, i, 3 ) = ai.X();
		    dydx ( sys, i, 4 ) = ai.Y();
		    dydx ( sys, i, 5 ) = ai.Z();
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
		const double ah[] = { 1.0 / 5.0, 0.3, 3.0 / 5.0, 1.0, 7.0 / 8.0 };
		const double b21 = 1.0 / 5.0;
		const double b3[] = { 3.0 / 40.0, 9.0 / 40.0 };
		const double b4[] = { 0.3, -0.9, 1.2 };
		const double b5[] = { -11.0 / 54.0, 2.5, -70.0 / 27.0, 35.0 / 27.0 };
		const double b6[] = { 1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0 };
		const double c1 = 37.0 / 378.0;
		const double c3 = 250.0 / 621.0;
		const double c4 = 125.0 / 594.0;
		const double c6 = 512.0 / 1771.0;
		const double ec[] = 
		     { 0.0,
		       /* the first value is the same as c1, above */
		       37.0 / 378.0 - 2825.0 / 27648.0,
		       0.0,
		       /* the first value is the same as c3, above */
		       250.0 / 621.0 - 18575.0 / 48384.0,
		       /* the first value is the same as c4, above */
		       125.0 / 594.0 - 13525.0 / 55296.0,
		       -277.00 / 14336.0,
		       /* the first value is the same as c6, above */
		       512.0 / 1771.0 - 0.25 };


			if(T >= Tend) { return T; }
//			double hh = T + this->h <= Tend ? this->h : Tend - T;

			// k1 step
			compute_acc(ens,sys,k1);       // at T
			double htmp = b21*h;
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + htmp* k1(sys, bod, 0);
			ytmp(sys, bod, 1) = ens.y (sys, bod) + htmp* k1(sys, bod, 1);
			ytmp(sys, bod, 2) = ens.z (sys, bod) + htmp* k1(sys, bod, 2);
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + htmp* k1(sys, bod, 3);
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + htmp* k1(sys, bod, 4);
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + htmp* k1(sys, bod, 5);
			}
			// k2 step
			compute_acc(ens,ytmp,sys,k2);  // at T + ah[0]*h
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + h* (b3[0]*k1(sys, bod, 0) + b3[1]*k2(sys, bod, 0));
			ytmp(sys, bod, 1) = ens.y (sys, bod) + h* (b3[0]*k1(sys, bod, 1) + b3[1]*k2(sys, bod, 1));
			ytmp(sys, bod, 2) = ens.z (sys, bod) + h* (b3[0]*k1(sys, bod, 2) + b3[1]*k2(sys, bod, 2));
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + h* (b3[0]*k1(sys, bod, 3) + b3[1]*k2(sys, bod, 3));
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + h* (b3[0]*k1(sys, bod, 4) + b3[1]*k2(sys, bod, 4));
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + h* (b3[0]*k1(sys, bod, 5) + b3[1]*k2(sys, bod, 5));
			}
			// k3 step
			compute_acc(ens,ytmp,sys,k3);  // at T + ah[1]*h
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + h* (b4[0]*k1(sys, bod, 0) + b4[1]*k2(sys, bod, 0) + b4[2]*k3(sys, bod, 0) );
			ytmp(sys, bod, 1) = ens.y (sys, bod) + h* (b4[0]*k1(sys, bod, 1) + b4[1]*k2(sys, bod, 1) + b4[2]*k3(sys, bod, 1) );
			ytmp(sys, bod, 2) = ens.z (sys, bod) + h* (b4[0]*k1(sys, bod, 2) + b4[1]*k2(sys, bod, 2) + b4[2]*k3(sys, bod, 2) );
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + h* (b4[0]*k1(sys, bod, 3) + b4[1]*k2(sys, bod, 3) + b4[2]*k3(sys, bod, 3) );
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + h* (b4[0]*k1(sys, bod, 4) + b4[1]*k2(sys, bod, 4) + b4[2]*k3(sys, bod, 4) );
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + h* (b4[0]*k1(sys, bod, 5) + b4[1]*k2(sys, bod, 5) + b4[2]*k3(sys, bod, 5) );
			}
			// k4 step
			compute_acc(ens,ytmp,sys,k4);// at T + ah[2]*h
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + h* (b5[0]*k1(sys, bod, 0) + b5[1]*k2(sys, bod, 0) + b5[2]*k3(sys, bod, 0) + b5[3]*k4(sys, bod, 0) );
			ytmp(sys, bod, 1) = ens.y (sys, bod) + h* (b5[0]*k1(sys, bod, 1) + b5[1]*k2(sys, bod, 1) + b5[2]*k3(sys, bod, 1) + b5[3]*k4(sys, bod, 1) );
			ytmp(sys, bod, 2) = ens.z (sys, bod) + h* (b5[0]*k1(sys, bod, 2) + b5[1]*k2(sys, bod, 2) + b5[2]*k3(sys, bod, 2) + b5[3]*k4(sys, bod, 2) );
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + h* (b5[0]*k1(sys, bod, 3) + b5[1]*k2(sys, bod, 3) + b5[2]*k3(sys, bod, 3) + b5[3]*k4(sys, bod, 3) );
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + h* (b5[0]*k1(sys, bod, 4) + b5[1]*k2(sys, bod, 4) + b5[2]*k3(sys, bod, 4) + b5[3]*k4(sys, bod, 4) );
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + h* (b5[0]*k1(sys, bod, 5) + b5[1]*k2(sys, bod, 5) + b5[2]*k3(sys, bod, 5) + b5[3]*k4(sys, bod, 5) );
			}
			// k5 step
			compute_acc(ens,ytmp,sys,k5);// at T + ah[3]*h
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			ytmp(sys, bod, 0) = ens.x (sys, bod) + h* (b6[0]*k1(sys, bod, 0) + b6[1]*k2(sys, bod, 0) + b6[2]*k3(sys, bod, 0) + b6[3]*k4(sys, bod, 0) + b6[4]*k5(sys, bod, 0));
			ytmp(sys, bod, 1) = ens.y (sys, bod) + h* (b6[0]*k1(sys, bod, 1) + b6[1]*k2(sys, bod, 1) + b6[2]*k3(sys, bod, 1) + b6[3]*k4(sys, bod, 1) + b6[4]*k5(sys, bod, 1));
			ytmp(sys, bod, 2) = ens.z (sys, bod) + h* (b6[0]*k1(sys, bod, 2) + b6[1]*k2(sys, bod, 2) + b6[2]*k3(sys, bod, 2) + b6[3]*k4(sys, bod, 2) + b6[4]*k5(sys, bod, 2));
			ytmp(sys, bod, 3) = ens.vx(sys, bod) + h* (b6[0]*k1(sys, bod, 3) + b6[1]*k2(sys, bod, 3) + b6[2]*k3(sys, bod, 3) + b6[3]*k4(sys, bod, 3) + b6[4]*k5(sys, bod, 3));
			ytmp(sys, bod, 4) = ens.vy(sys, bod) + h* (b6[0]*k1(sys, bod, 4) + b6[1]*k2(sys, bod, 4) + b6[2]*k3(sys, bod, 4) + b6[3]*k4(sys, bod, 4) + b6[4]*k5(sys, bod, 4));
			ytmp(sys, bod, 5) = ens.vz(sys, bod) + h* (b6[0]*k1(sys, bod, 5) + b6[1]*k2(sys, bod, 5) + b6[2]*k3(sys, bod, 5) + b6[3]*k4(sys, bod, 5) + b6[4]*k5(sys, bod, 5));
			}
			// k6 step & final sum
			compute_acc(ens,ytmp,sys,k6);// at T + ah[4]*h
			for(int bod = 0; bod != ens.nbod(); bod++)
			{
			double di;
			di = c1*k1(sys, bod, 0) + c3*k3(sys, bod, 0) + c4*k4(sys, bod, 0) + c6*k6(sys, bod, 0);
			ens.x(sys,bod) += h * di;
			//  yerr(sys, bod, 0) = h* ( ec[1]*k1(sys, bod, 0) + ec[3]*k3(sys, bod, 0)  ec[4] * k4(sys, bod, 0) + ec[5] * k5(sys, bod, 0) + ec[6] * k6(sys, bod, 0) );
			di = c1*k1(sys, bod, 1) + c3*k3(sys, bod, 1) + c4*k4(sys, bod, 1) + c6*k6(sys, bod, 1);
			ens.y(sys,bod) += h * di;
			//  yerr(sys, bod, 1) = h* ( ec[1]*k1(sys, bod, 1) + ec[3]*k3(sys, bod, 1)  ec[4] * k4(sys, bod, 1) + ec[5] * k5(sys, bod, 1) + ec[6] * k6(sys, bod, 1) );
			di = c1*k1(sys, bod, 2) + c3*k3(sys, bod, 2) + c4*k4(sys, bod, 2) + c6*k6(sys, bod, 2);
			ens.z(sys,bod) += h * di;
			//  yerr(sys, bod, 2) = h* ( ec[1]*k1(sys, bod, 2) + ec[3]*k3(sys, bod, 2)  ec[4] * k4(sys, bod, 2) + ec[5] * k5(sys, bod, 2) + ec[6] * k6(sys, bod, 2) );
			di = c1*k1(sys, bod, 3) + c3*k3(sys, bod, 3) + c4*k4(sys, bod, 3) + c6*k6(sys, bod, 3);
			ens.vx(sys,bod) += h * di;
			//  yerr(sys, bod, 3) = h* ( ec[1]*k1(sys, bod, 3) + ec[3]*k3(sys, bod, 3)  ec[4] * k4(sys, bod, 3) + ec[5] * k5(sys, bod, 3) + ec[6] * k6(sys, bod, 3) );
			di = c1*k1(sys, bod, 4) + c3*k3(sys, bod, 4) + c4*k4(sys, bod, 4) + c6*k6(sys, bod, 4);
			ens.vy(sys,bod) += h * di;
			//  yerr(sys, bod, 4) = h* ( ec[1]*k1(sys, bod, 4) + ec[3]*k3(sys, bod, 4)  ec[4] * k4(sys, bod, 4) + ec[5] * k5(sys, bod, 4) + ec[6] * k6(sys, bod, 4) );
			di = c1*k1(sys, bod, 5) + c3*k3(sys, bod, 5) + c4*k4(sys, bod, 5) + c6*k6(sys, bod, 5);
			ens.vz(sys,bod) += h * di;
			//  yerr(sys, bod, 5) = h* ( ec[1]*k1(sys, bod, 5) + ec[3]*k3(sys, bod, 5)  ec[4] * k4(sys, bod, 5) + ec[5] * k5(sys, bod, 5) + ec[6] * k6(sys, bod, 5) );
			}

			return T + h;
		}

	};

	//! CPU state and interface
	cuxDeviceAutoPtr<double, 3> k1, k2, k3, k4, k5, k6;
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
		k1.realloc(ens.nsys(), ens.nbod(), 6);
		k2.realloc(ens.nsys(), ens.nbod(), 6);
		k3.realloc(ens.nsys(), ens.nbod(), 6);
		k4.realloc(ens.nsys(), ens.nbod(), 6);
		k5.realloc(ens.nsys(), ens.nbod(), 6);
		k6.realloc(ens.nsys(), ens.nbod(), 6);
//		y0.realloc(ens.nsys(), ens.nbod(), 6);
		ytmp.realloc(ens.nsys(), ens.nbod(), 6);
//		yerr.realloc(ens.nsys(), ens.nbod(), 6);
		gpu_obj.k1= k1;
		gpu_obj.k2= k2;
		gpu_obj.k3= k3;
		gpu_obj.k4= k4;
		gpu_obj.k5= k5;
		gpu_obj.k6= k6;
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
	prop_rkck(const config &cfg)
	{
		if(!cfg.count("time step")) ERROR("Integrator gpu_rkck requires a timestep ('time step' keyword in the config file).");
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
 * with the propagator rkck and the stopper stop_on_ejection
 * 
 * @param[in] cfg contains configuration data for gpu_generic_integrator
 */
extern "C" integrator *create_gpu_rkck(const config &cfg)
{
	return new gpu_generic_integrator<stop_on_ejection, prop_rkck>(cfg);
}

} // end namespace swarm
