/*************************************************************************
 * Copyright (C) 2011 by Eric B. Ford and the Swarm-NG Development Team  *
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

// #define N_LAG 5.0 // integer n, for recommended Laguerre method
#define SIGN(a) ((a) < 0 ? -1 : 1)


namespace swarm {
namespace hp {

class mvs: public integrator {
	typedef integrator base;
	private:
	double _time_step;
	static const int _N_LAG = 5.0;

	public:
	mvs(const config& cfg): base(cfg),_time_step(0.001) {
		if(!cfg.count("time step")) ERROR("Integrator gpu_mvs requires a timestep ('time step' keyword in the config file).");
		_time_step = atof(cfg.at("time step").c_str());
	}

	virtual void launch_integrator() {
		launch_templatized_integrator(this);
	}



////////////////////////////////////////////////////////////////
// solving differential Kepler's equation
// in universal variable x
// using Laguerre method as outlined by Prusing+C eqn 2.43
// code adapted from Alice Quillen's Qymsym code 
// see http://astro.pas.rochester.edu/~aquillen/qymsym/
////////////////////////////////////////////////////////////////
#define MINR 1.0e-5 // minimum radius
#define MINDENOM 1e-8  // mininum denominator

__device__ double solvex(double r0dotv0, double alpha,
                double M1, double r0, double dt)
{
   double smu = sqrt(M1);
   double foo = 1.0 - r0*alpha;
   double sig0 = r0dotv0/smu;
   double x = M1*dt*dt/r0; // initial guess could be improved 
// better initial guess depends on rperi which would have to be passed

   double u=1.0;
   for(int i=0;i<7;i++){  // 7 iterations is probably overkill
			// as it always converges faster than this
     double x2,x3,alx2,Cp,Sp,F,dF,ddF,z;
     x2 = x*x;
     x3 = x2*x;
     alx2 = alpha*x2;
     Cp = C_prussing(alx2);
     Sp = S_prussing(alx2);
     F = sig0*x2*Cp + foo*x3*Sp + r0*x - smu*dt; // eqn 2.41 PC
     dF = sig0*x*(1.0 - alx2*Sp)  + foo*x2*Cp + r0; // eqn 2.42 PC
     ddF = sig0*(1.0-alx2*Cp) + foo*x*(1.0 - alx2*Sp);
     z = fabs((_N_LAG - 1.0)*((_N_LAG - 1.0)*dF*dF - _N_LAG*F*ddF));
     z = sqrt(z);
     double denom = (dF + SIGN(dF)*z); 
     if (denom ==0.0) denom = MINDENOM;
     u = _N_LAG*F/denom; // equation 2.43 PC
     x -= u;
   }
//   if (isnan(x)) printf("solvex: is nan\n");
   return x;
}

// functions needed for kepstep
// code adapted from Alice Quillen's Qymsym code 
// see http://astro.pas.rochester.edu/~aquillen/qymsym/
__device__ double C_prussing(double y) // equation 2.40a Prussing + Conway
{
  if (fabs(y)<1e-4) return 1.0/2.0*(1.0 - y/12.0*(1.0 - y/30.0*(1.0 - y/56.0)));
  double u = sqrt(fabs(y));
  if (y>0.0) return (1.0- cos(u))/ y;
  else       return (cosh(u)-1.0)/-y;
}

__device__ double S_prussing(double y) // equation 2.40b Prussing +Conway
{
  if (fabs(y)<1e-4) return 1.0/6.0*(1.0 - y/20.0*(1.0 - y/42.0*(1.0 - y/72.0)));
  double u = sqrt(fabs(y));
  double u3 = u*u*u;
  if (y>0.0) return (u -  sin(u))/u3;
  else       return (sinh(u) - u)/u3;
}

///////////////////////////////////////////////////////////////
// advance a particle using f,g functions and universal variables
// for differental kepler's equation
//  has an error catch for r0=0 so can be run with central body
// Based on equations by Prussing, J. E. and Conway, B. A. 
// Orbital Mechanics 1993, Oxford University Press, NY,NY  chap 2 
// npos,nvel are new positions and velocity
// pos, vel are old ones
// code adapted from Alice Quillen's Qymsym code 
// see http://astro.pas.rochester.edu/~aquillen/qymsym/
///////////////////////////////////////////////////////////////
//__device__ void kepstep(double4 pos, double4 vel, double4* npos, double4* nvel, double deltaTime, double GM)
__device__ void drift_kepler(double& x, double& y, double& z, double& vx, double& vy, double& vz, const double GM, const double deltaTime)
{
   double x_old = x, y_old = y, z_old = z, vx_old = vx, vy_old = vy, vz_old = vz;
   // WARNING: Using softened potential
   double r0 = sqrt(x*x + y*y + z*z + MINR*MINR); // current radius
   double v2 = (vx*vx + vy*vy + vz*vz);  // current velocity
   double r0dotv0 = (x*vx + y*vy + z*vz);
   double alpha = (2.0/r0 - v2/GM);  // inverse of semi-major eqn 2.134 MD
// here alpha=1/a and can be negative
   double x_p = solvex(r0dotv0, alpha, GM, r0, deltaTime); // solve universal kepler eqn

   double smu = sqrt(GM);  // TODO: probably could cache sqrt(GM)
   double foo = 1.0 - r0*alpha;
   double sig0 = r0dotv0/smu;
   double x2 = x_p*x_p;
   double x3 = x2*x_p;
   double alx2 = alpha*x2;
   double Cp = C_prussing(alx2);
   double Sp = S_prussing(alx2);
   double r = sig0*x_p*(1.0 - alx2*Sp)  + foo*x2*Cp + r0; // eqn 2.42  PC
   if (r < MINR) r=MINR;
// if dt == 0 then f=dgdt=1 and g=dfdt=0
// f,g functions equation 2.38a  PC
   double f_p= 1.0 - (x2/r0)*Cp;
   double g_p= deltaTime - (x3/smu)*Sp;
// dfdt,dgdt function equation 2.38b PC
   double dfdt;
   double dgdt = 1.0 - (x2/r)*Cp;
   if (fabs(g_p) > MINDENOM)
      // conservation of angular momentum means that f dfdt - g dfdt =1
      dfdt = (f_p*dgdt - 1.0)/g_p;
   else
      // dfdt,dgdt function equation 2.38b PC
      dfdt = x_p*smu/(r*r0)*(alx2*Sp - 1.0);
  
   x = f_p*x_old + g_p*vx_old;     // eqn 2.65 M+D
   y = f_p*y_old + g_p*vy_old;
   z = f_p*z_old + g_p*vz_old; 
   vx = dfdt*x_old + dgdt*vx_old; //eqn 2.70 M+D
   vy = dfdt*y_old + dgdt*vy_old;
   vz = dfdt*z_old + dgdt*vz_old;

}

	/*! Integrator Kernel to be run on GPU
	 *  
	 *
	 */
	 template<class T >
	__device__ void kernel(T a)  {

		const int nbod = T::n;

		/////////////// FETCH LOCAL VARIABLES ///////////////////

		int thr = thread_in_system();

		/* WARNING: This is fine for now, but I wonder if this should be removed in the future.
		   I think it's here to prevent readingn outside the bounds of _gpu_ens.
		   But I worry that if we're usingn the threads in different ways, are we sure we'd never want to use this thread for something?  */
		if(sysid() >= _gpu_ens->nsys()) { return; }

		ensemble::systemref sys ( (*_gpu_ens)[sysid()] );

		// Body/Component Grid
		// Body number
		int b = thr / 3 ;
		int bb = b+1; // for excluding central body
		// Component number
		int c = thr % 3 ;
//		bool body_component_grid = b < nbod;
		bool body_component_grid_no_sun = bb < nbod;
//		bool body_grid = thr < nbod;

		// i,j pairs Grid
		int ij = thr;

		// shared memory allocation
		extern __shared__ char shared_mem[];
		char*  system_shmem =( shared_mem + sysid_in_block() * integrator::shmem_per_system(nbod) );

		double t_start = sys.time(), t = t_start;
		double t_end = min(t_start + _destination_time,sys.time_end());

		// local information per component per body
		double pos = 0, vel = 0 , acc = 0, jerk = 0;
		if( body_component_grid_no_sun )
			pos = sys[bb].p(c), vel = sys[bb].v(c);

		////////// INTEGRATION //////////////////////

		// Calculate acceleration and jerk
		Gravitation<nbod> calcForces(sys,system_shmem);
		calcForces.calc_accel_no_sun(ij,bb,c,acc,jerk);

		unsigned int iter=0;
		while(t < t_end)
		{
		   double hby2 = 0.5*min(_time_step, t_end - t);
//		   double pos_old = pos, vel_old = vel, acc_old = acc,jerk_old = jerk;	

		   // Kick
		   if( body_component_grid_no_sun )
		      {
//		      pos = pos +  hby2*(vel+(hby2*0.5)*(acc+(hby2/3.)*jerk));
//		      sys[bb].p(c) = pos; 
		      vel = vel +  hby2*(acc+(hby2*0.5)*jerk);
		      sys[bb].v(c) = vel;
		      }
		   __syncthreads();

		   // Drift
		   int bbb = thr+1;
		   if(bbb < nbod)  // Central body does not drift
		   {
		   // TODO: For now using heliocentric drift; improve
		   // TODO: Use algorithm that works; currently in non-inertial frame, but not including ficticous forces to compensate
		   double cx = sys[0].p(0), cy = sys[0].p(1), cz = sys[0].p(2);
		   double cvx = sys[0].v(0), cvy = sys[0].v(1), cvz = sys[0].v(2);
		   double x = sys[bbb].p(0)-cx, y = sys[bbb].p(1)-cy, z = sys[bbb].p(2)-cz;
		   double vx = sys[bbb].v(0)-cvx, vy = sys[bbb].v(1)-cvy, vz = sys[bbb].v(2)-cvz;
		   // TODO: could probably cache sqrt(mass) across steps
		   double GM = sys[0].mass()+sys[bbb].mass();
                   drift_kepler(x,y,z,vx,vy,vz,GM,2.0*hby2);
		   sys[bbb].p(0) = cx + x; // + 2.0*hby2*cvx;
		   sys[bbb].p(1) = cy + y; // + 2.0*hby2*cvy;
		   sys[bbb].p(2) = cz + z; // + 2.0*hby2*cvz;
		   sys[bbb].v(0) = cvx + vx;
		   sys[bbb].v(1) = cvy + vy;
		   sys[bbb].v(2) = cvz + vz;
		   }
		   __syncthreads();	   		   

		   // Kick
		   calcForces.calc_accel_no_sun(ij,bb,c,acc,jerk);

		   if( body_component_grid_no_sun )
		     {
//		     pos = pos +  hby2*(vel+(hby2*0.5)*(acc+(hby2/3.)*jerk));
//		     sys[bb].p(c) = pos;
		     vel = vel +  hby2*(acc+(hby2*0.5)*jerk);
		     sys[bb].v(c) = vel;
		     }
		   __syncthreads(); // in case will output from sys


		   t += 2.*hby2;
		   ++iter;
		   if((thr == 0)  && (log::needs_output(*_gpu_ens, t, sysid())) )
		      {
		      sys.set_time(t);
		      log::output_system(*_gpu_log, *_gpu_ens, t, sysid());
		      }

		   if(iter>=_max_itterations_per_kernel_call) break;
		}

		if(thr == 0) 
		   sys.set_time(t);

	}

};

/*!
 * \brief Factory to create mvs gpu integrator from config class
 *
 * @param[in] cfg configuration class
 *
 * @return        pointer to integrator cast to integrator*
 */
/*
extern "C" integrator *create_hp_mvs_fixed(const config &cfg)
{
	return new mvs< FixedTimeStep> (cfg);
}

extern "C" integrator *create_hp_mvs_adaptive(const config &cfg)
{
	return new mvs< AdaptiveTimeStep> (cfg);
}
*/
extern "C" integrator *create_hp_mvs(const config &cfg)
{
	return new mvs(cfg);
}

}
}

