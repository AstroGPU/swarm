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
                double sqrtM1, double r0, double dt)
{
//   double smu = sqrt(M1);
   double smu = sqrtM1;
   double foo = 1.0 - r0*alpha;
   double sig0 = r0dotv0/smu;
   double x = sqrtM1*sqrtM1*dt*dt/r0; // initial guess could be improved 
// better initial guess depends on rperi which would have to be passed

   double u=1.0;
   // TODO: Could allow to exit loop early
   for(int i=0;i<7;i++){  // 7 iterations is probably overkill
			// as it always converges faster than this
     double x2,x3,alx2,Cp,Sp,F,dF,ddF,z;
     x2 = x*x;
     x3 = x2*x;
     alx2 = alpha*x2;
     Cp = C_prussing(alx2);
     Sp = S_prussing(alx2);
//   TODO: Test if this is faster/accurate     
//   SC_prussing(alx2,Sp,Cp);
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
// TODO: could merge into one function that uses sincos
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

__device__ void SC_prussing(double y, double& S, double &C) // equation 2.40a Prussing + Conway
{
  if (fabs(y)<1e-4) 
  {
     S = 1.0/6.0*(1.0 - y/20.0*(1.0 - y/42.0*(1.0 - y/72.0)));
     C = 1.0/2.0*(1.0 - y/12.0*(1.0 - y/30.0*(1.0 - y/56.0)));
     return;
  }
  double u = sqrt(fabs(y));
  double u3 = u*u*u;
  if (y>0.0) 
     {
     sincos(u,&S,&C);  // TODO: Need to verify called correctly
     S = (u -  S)/u3;
     C = (1.0- C)/ y;
     }
  else
     {
       S = (sinh(u) - u)/u3;
       C = (cosh(u)-1.0)/-y;
     }
  return;
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
__device__ void drift_kepler(double& x_old, double& y_old, double& z_old, double& vx_old, double& vy_old, double& vz_old, const double sqrtGM, const double deltaTime)
{
   double x = x_old, y = y_old, z = z_old, vx = vx_old, vy = vy_old, vz = vz_old;
   // WARNING: Using softened potential
   double r0 = sqrt(x*x + y*y + z*z + MINR*MINR); // current radius
   double v2 = (vx*vx + vy*vy + vz*vz);  // current velocity
   double r0dotv0 = (x*vx + y*vy + z*vz);
   double GM = sqrtGM*sqrtGM;
   double alpha = (2.0/r0 - v2/GM);  // inverse of semi-major eqn 2.134 MD
// here alpha=1/a and can be negative
   double x_p = solvex(r0dotv0, alpha, sqrtGM, r0, deltaTime); // solve universal kepler eqn

//   double smu = sqrt(GM);  // before we cached sqrt(GM)
   double smu = sqrtGM; 
   double foo = 1.0 - r0*alpha;
   double sig0 = r0dotv0/smu;
   double x2 = x_p*x_p;
   double x3 = x2*x_p;
   double alx2 = alpha*x2;
   double Cp = C_prussing(alx2);
   double Sp = S_prussing(alx2);
//   TODO: Test if this is faster/accurate     
//   double Cp, Sp;
//   SC_prussing(alx2,Sp,Cp);
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
   vx = dfdt*x_old + dgdt*vx_old;  // eqn 2.70 M+D
   vy = dfdt*y_old + dgdt*vy_old;
   vz = dfdt*z_old + dgdt*vz_old;

   // Replace values 
    x_old =  x;  y_old =  y;  z_old =  z;
   vx_old = vx; vy_old = vy; vz_old = vz;
}

	/*! Integrator Kernel to be run on GPU
	 *  
	 * TODO: Need to deal with energy conservation if input not in COM frame
	 */
	 template<class T >
	__device__ void kernel(T a)  {

		const int nbod = T::n;
		const bool allow_rewind = false;

		/////////////// FETCH LOCAL VARIABLES ///////////////////

		int thr = thread_in_system();

		/* WARNING: This is fine for now, but I wonder if this should be removed in the future.
		   I think it's here to prevent readingn outside the bounds of _gpu_ens.
		   But I worry that if we're usingn the threads in different ways, are we sure we'd never want to use this thread for something?  */
		if(sysid() >= _gpu_ens->nsys()) { return; }

		ensemble::systemref sys ( (*_gpu_ens)[sysid()] );

		// Body/Component Grid
		// Body number
		int b = thr / 3 ;  // index for parts w/ 1 thread per body
		int bb = b+1;      // index for parts w/ 1 thread per body excluding sun/central body
		// Component number
		int c = thr % 3 ;  
		bool body_component_grid = b < nbod;          // if needed for parts w/ 1 thread per body per component including sun/central body
		bool body_component_grid_no_sun = bb < nbod;  // if needed for parts w/ 1 thread per body per component excluding sun/central body
//		bool body_grid = thr < nbod;                  // if needed for parts w/ 1 thread per body including sun/central body

		// i,j pairs Grid
		// TODO: Be more clever about calculating accelerations
		//       Either avoid pairs with sun or specialize Gravitation?
		int ij = thr;      // index for parts w/ 1 thread per body pair


		// shared memory allocation
		extern __shared__ char shared_mem[];
		char*  system_shmem =( shared_mem + sysid_in_block() * integrator::shmem_per_system(nbod) );

		double t_start = sys.time(), t = t_start;
		double t_end = min(t_start + _destination_time,sys.time_end());

		// local information per component per body
		double pos_old, vel_old, acc_old, jerk_old; // needed if allowing rewindss
		double acc = 0., jerk = 0.;
		double sqrtGM = sqrt(sys[0].mass()); // TODO: Could parallelize. Worth it?

		// Shift into funky coordinate system (see A. Quillen's qymsym's tobary)
		if( (b==0) || body_component_grid_no_sun )
		   {
		   double sump = 0., sumv = 0., mtot = 0.;
		   for(int j=0;j<nbod;++j)   // TODO: Could parallelize. Worth it?
		      {
		      const double mj = sys[j].mass();
		      mtot += mj;
		      sump += mj*sys[j].p(c);
		      sumv += mj*sys[j].v(c);
		      }
		   if(b==0) // For sun only
		      {
		      sys[b].v(c) = sumv/mtot;
		      sys[b].p(c) = sump/mtot;
		      }
		   if( body_component_grid_no_sun ) // For all bodies except sun
		      {
		      sys[bb].v(c) -= sumv/mtot;
		      sys[bb].p(c) -= sys[0].p(c);
  		      }
		   }
		   __syncthreads();		

		////////// INTEGRATION //////////////////////
		// Calculate acceleration and jerk
		Gravitation<nbod> calcForces(sys,system_shmem);
                // precompute acc and jerk before enter loop, since will cache acc and jerk across loop itterations
		calcForces.calc_accel_no_sun(ij,bb,c,acc,jerk);

		unsigned int iter=0;  // Make sure don't get stuck in infinite loop
		while(t < t_end)      // Only enter loop if need to integrate
		{
		   double hby2 = 0.5*min(_time_step, t_end - t);
		   if(allow_rewind)   // Could be useful if later reject step
		     { 
		     if( body_component_grid )
			{ 
			pos_old = sys[b].p(c); vel_old = sys[b].v(c); 
			acc_old = acc;  jerk_old = jerk;
			}
		     }

		   // Drift Step (center-of-mass motion)
		   if( body_component_grid_no_sun )
		      {
		      double mv = 0.;
		      // TODO: In principle could parellalize.  Worth it?
		      for(int j=1;j<nbod;++j)
		      	 mv += sys[j].mass()*sys[j].v(c);
		      sys[bb].p(c) += mv*hby2/sys[0].mass();
		      }
		   __syncthreads();

		   // Kick Step (planet-planet interactions)
		   {
		   // TODO: Test that this call can be removed
		   // WARNING: If make changes, check that it's ok to not recompute
		   calcForces.calc_accel_no_sun(ij,bb,c,acc,jerk);
		   if( body_component_grid_no_sun )
		      {
		      sys[bb].v(c) +=  hby2*(acc+hby2*0.5*jerk);
		      }
		   }
		   __syncthreads();
  
		   // Kepler Drift Step (Keplerian orbit about sun/central body)
		   int bbb = thr+1;
		   if(bbb < nbod)  // Central body does not do Kepler drift
		   {
		      drift_kepler(sys[bbb].p(0),sys[bbb].p(1),sys[bbb].p(2),sys[bbb].v(0),sys[bbb].v(1),sys[bbb].v(2),sqrtGM,2.0*hby2);
		   }
		   __syncthreads();	   		   
		   
		   /* TODO: Eventually check for close encounters and 
		            if necessary undo, perform direct n-body, merge and resume
	  	            Or maybe only in separate integrator? */
		   bool need_to_rewind = false;
		   if( allow_rewind && need_to_rewind )
		     {
			sys[b].p(c) = pos_old; sys[b].v(c) = vel_old; 
			acc = acc_old;  jerk = jerk_old;
		     	++iter;
		   	if(iter>=_max_itterations_per_kernel_call) break;
			continue;
                     }

		   // Kick Step (planet-planet interactions)
		   {
		   calcForces.calc_accel_no_sun(ij,bb,c,acc,jerk);
		   if( body_component_grid_no_sun )
		      {
		      sys[bb].v(c) +=  hby2*(acc+hby2*0.5*jerk);
		      }
		   }
		   __syncthreads();

		   // Drift Step (center-of-mass motion)
		   if( body_component_grid_no_sun )
		      {
		      double mv = 0.;
		      // TODO: In principle could parellalize. Worth it?
		      for(int j=1;j<nbod;++j)
		      	 mv += sys[j].mass()*sys[j].v(c);
		      sys[bb].p(c) += mv*hby2/sys[0].mass();
		      }
		   __syncthreads();

		   // WARNING: Need to think about correct order of time updates, if add time dependnt forces
		   t += 2.*hby2;

		   ++iter;
		   if( log::needs_output(*_gpu_ens, t, sysid()) )
		      {
		      // Save working coordinates
		      double pos_tmp, vel_tmp;
		      if(body_component_grid )
			{ pos_tmp = sys[b].p(c); vel_tmp = sys[b].v(c); }

		      // Shift back from funky coordinate system (see A. Quillen's qymsym's tobary)
		      if( (b==0) || body_component_grid_no_sun )
		         {
		   	 const double m0 = sys[0].mass();
		   	 double sump = 0., sumv = 0., mtot = m0;
		   	 for(int j=1;j<nbod;++j)   // TODO: Could parallelize;  Worth it?
		      	    {
		      	    const double mj = sys[j].mass();
		      	    mtot += mj;
		      	    sump += mj*sys[j].p(c);
		      	    sumv += mj*sys[j].v(c);
		      	    }
		   	 if(b==0) // For sun only
		      	    {
		      	    sys[b].p(c) -= sump/mtot;
		      	    sys[b].v(c) -= sumv/m0;
		      	    }
		   	 if( body_component_grid_no_sun ) // For all bodies except sun
		      	    {
		      	    sys[bb].p(c) += sys[0].p(c) - sump/mtot;
		      	    sys[bb].v(c) += sys[0].v(c);
  		      	    }
		   	 }
		      __syncthreads();		      
		      if(thr == 0)
		         {
		         sys.set_time(t);
		         log::output_system(*_gpu_log, *_gpu_ens, t, sysid());
		         }
		      __syncthreads();
		      // Restore working coordinates
		      if(body_component_grid )
			{ sys[b].p(c) = pos_tmp; sys[b].v(c) = vel_tmp; }
		      __syncthreads();
		      }

		   if(iter>=_max_itterations_per_kernel_call) break;
		}

		// Shift back from funky coordinate system (see A. Quillen's qymsym's tobary)
		if( (b==0) || body_component_grid_no_sun )
		   {
		   const double m0 = sys[0].mass();
		   double sump = 0., sumv = 0., mtot = m0;
		   for(int j=1;j<nbod;++j)   // TODO: Could parallelize. Worth it?
		      {
		      const double mj = sys[j].mass();
		      mtot += mj;
		      sump += mj*sys[j].p(c);
		      sumv += mj*sys[j].v(c);
		      }
		   if(b==0) // For sun only
		      {
		      sys[b].p(c) -= sump/mtot;
		      sys[b].v(c) -= sumv/m0;
		      }
		   if( body_component_grid_no_sun ) // For all bodies except sun
		      {
		      sys[bb].p(c) += sys[0].p(c) - sump/mtot;
		      sys[bb].v(c) += sys[0].v(c);
  		      }
		   }
		   __syncthreads();

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

