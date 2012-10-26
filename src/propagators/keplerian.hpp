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

/*! \file keplerian.hpp
 *   \brief Defines a solver for differential Kepler's equation in universal variable x. 
 *
 */


////////////////////////////////////////////////////////////////
// solving differential Kepler's equation
// in universal variable x
// using Laguerre method as outlined by Prusing+C eqn 2.43
// code adapted from Alice Quillen's Qymsym code 
// see http://astro.pas.rochester.edu/~aquillen/qymsym/
////////////////////////////////////////////////////////////////

#ifndef H_KEPLERIAN
#define H_KEPLERIAN

//#define MINR 1.0e-5 // minimum radius
#define MINR_IN_1EM8 0 // minimum radius
#define MINDENOM 1e-8  // mininum denominator
#define SIGN(a) ((a) < 0 ? -1 : 1)

// functions needed for kepstep
// code adapted from Alice Quillen's Qymsym code 
// see http://astro.pas.rochester.edu/~aquillen/qymsym/
GPUAPI double C_prussing(double y) // equation 2.40a Prussing + Conway
{
  if (fabs(y)<1e-4) return 1.0/2.0*(1.0 - y/12.0*(1.0 - y/30.0*(1.0 - y/56.0)));
  double u = sqrt(fabs(y));
  if (y>0.0) return (1.0- cos(u))/ y;
  else       return (cosh(u)-1.0)/-y;
}

GPUAPI double S_prussing(double y) // equation 2.40b Prussing +Conway
{
  if (fabs(y)<1e-4) return 1.0/6.0*(1.0 - y/20.0*(1.0 - y/42.0*(1.0 - y/72.0)));
  double u = sqrt(fabs(y));
  double u3 = u*u*u;
  if (y>0.0) return (u -  sin(u))/u3;
  else       return (sinh(u) - u)/u3;
}

GPUAPI void SC_prussing(double y, double& S, double &C) // equation 2.40a Prussing + Conway
{
  if (fabs(y)<1e-4) 
     {
       S = (1.0/6.0)*(1.0 - y*(1.0/20.0)*(1.0 - y*(1.0/42.0)*(1.0 - y*(1.0/72.0))));
       C = (1.0/2.0)*(1.0 - y*(1.0/12.0)*(1.0 - y*(1.0/30.0)*(1.0 - y*(1.0/56.0))));
     }
  else
     {
       float absy = y;
       absy = fabs(y);
     double u = sqrt(absy);
     double u3 = u*u*u;
     if (y>0.0) 
        {
     	sincos(u,&S,&C); 
     	S = (u -  S)/u3;
     	C = (1.0- C)/ y;
     	}
     else
	{
     	S = (sinh(u) - u)/u3;
     	C = (cosh(u)-1.0)/-y;
     	}
     }
  return;
}

__device__ void SC_prussing_fast(double y, double& S, double &C) // equation 2.40a Prussing + Conway
{
  if (y*y<1e-8) 
     {
       S = (1.0/6.0)*(1.0 - y*(1.0/20.0)*(1.0 - y*(1.0/42.0)*(1.0 - y*(1.0/72.0))));
       C = (1.0/2.0)*(1.0 - y*(1.0/12.0)*(1.0 - y*(1.0/30.0)*(1.0 - y*(1.0/56.0))));
     }
  else
     {
     float absy = y;
     absy = fabs(y);
     float u = sqrtf(absy);
     float u3 = u*u*u;
     if (y>0.0) 
        {
	  float Sf, Cf;
	  sincosf(u,&Sf,&Cf); 
	  S = (u -  Sf)/u3;
	  C = (1.0- Cf)/ y;
     	}
     else
	{
	  S = (sinhf(u) - u)/u3;
	  C = (coshf(u)-1.0)/-y;
     	}
     }
  return;
}


GPUAPI double solvex(double r0dotv0, double alpha,
                double sqrtM1, double r0, double dt)
{
   const double _N_LAG = 5.0; // integer n, for recommended Laguerre method
//   double smu = sqrt(M1);
   double smu = sqrtM1;
   double foo = 1.0 - r0*alpha;
   double sig0 = r0dotv0/smu;
   double x = sqrtM1*sqrtM1*dt*dt/r0; // initial guess could be improved 
// better initial guess depends on rperi which would have to be passed

   double u=1.0;
   for(int i=0;(i<7)&&!((i>2)&&(x+u==x));i++){  // 7 iterations is probably overkill
			  // as it always converges faster than this
     double x2,x3,alx2,Cp,Sp,F,dF,ddF,z;
     x2 = x*x;
     x3 = x2*x;
     alx2 = alpha*x2;
     // Cp = C_prussing(alx2);
     // Sp = S_prussing(alx2);
#if 0
     if(i==0)
       SC_prussing_fast(alx2,Sp,Cp);  // optimization
     else
#endif
     SC_prussing(alx2,Sp,Cp);  // optimization
     F = sig0*x2*Cp + foo*x3*Sp + r0*x - smu*dt; // eqn 2.41 PC
     dF = sig0*x*(1.0 - alx2*Sp)  + foo*x2*Cp + r0; // eqn 2.42 PC
     ddF = sig0*(1.0-alx2*Cp) + foo*x*(1.0 - alx2*Sp);
     z = fabs((_N_LAG - 1.0)*((_N_LAG - 1.0)*dF*dF - _N_LAG*F*ddF));
     z = sqrt(z);
     double denom = (dF + SIGN(dF)*z);  // faster than copysign
     if (denom ==0.0) denom = MINDENOM;
     u = _N_LAG*F/denom; // equation 2.43 PC
     x -= u;
//     if( (i>=2) && (x+u==x) ) break;      // optimization to allow to exit loop early
   }
//   if (isnan(x)) printf("solvex: is nan\n");
   return x;
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
//GPUAPI void kepstep(double4 pos, double4 vel, double4* npos, double4* nvel, double deltaTime, double GM)
GPUAPI void drift_kepler(double& x_old, double& y_old, double& z_old, double& vx_old, double& vy_old, double& vz_old, const double sqrtGM, const double deltaTime)
{
   double x = x_old, y = y_old, z = z_old, vx = vx_old, vy = vy_old, vz = vz_old;
#if (MINR_IN_1EM8>0)
   // WARNING: Using softened potential
   double r0 = sqrt(x*x + y*y + z*z + MINR_IN_1EM8*MINR_IN_1EM8*1.e-16); // current radius
#else
   double r0 = sqrt(x*x + y*y + z*z ); // current radius
#endif
//   double r0 = sqrt(x*x + y*y + z*z ); // current radius
   double v2 = (vx*vx + vy*vy + vz*vz);  // current velocity
   double r0dotv0 = (x*vx + y*vy + z*vz);
   double GM = sqrtGM*sqrtGM;
   double alpha = (2.0/r0 - v2/GM);  // inverse of semi-major eqn 2.134 MD
// here alpha=1/a and can be negative
   double x_p = solvex(r0dotv0, alpha, sqrtGM, r0, deltaTime); // solve universal kepler eqn

//   double smu = sqrt(GM);  // from before we cached sqrt(GM)
   double smu = sqrtGM; 
   double foo = 1.0 - r0*alpha;
   double sig0 = r0dotv0/smu;
   double x2 = x_p*x_p;
   double x3 = x2*x_p;
   double alx2 = alpha*x2;
//   double Cp = C_prussing(alx2);
//   double Sp = S_prussing(alx2);
   double Cp, Sp;
   SC_prussing(alx2,Sp,Cp);  // optimization
   double r = sig0*x_p*(1.0 - alx2*Sp)  + foo*x2*Cp + r0; // eqn 2.42  PC

#if (MINR_IN_1EM8>0)
// WARNING: Using softened potentil
   if (r < MINR_IN_1EM8*1.e-8) r=MINR_IN_1EM8*1.e-8;
#else
   if(r<0.) r = 0.;
#endif

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

#endif
