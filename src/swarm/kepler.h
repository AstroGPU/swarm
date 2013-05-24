/*************************************************************************
 * Copyright (C) 2009-2010 by Eric Ford & the Swarm-NG Development Team  *
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


/*! \file kepler.h
 *  \brief converts between Cartesian & Keplerian coordinates
 *
 *  Based on code from John Chambers' Mercury code 
 *  \todo put proper refernece to publication or URL.
 */

#pragma once

double improve_mean_to_eccentric_annomaly_guess(const double e, const double M, const double x);
double mean_to_eccentric_annomaly(const double e,  double M);
void calc_cartesian_for_ellipse(double& x,double& y, double & z, double &vx, double &vy, double &vz, const double a, const double e, const double i, const double O, const double w, const double M, const double GM);
void calc_keplerian_for_cartesian( double& a,  double& e,  double& i,  double& O,  double& w,  double& M, const double x,const double y, const double z, const double vx, const double vy, const double vz, const double GM);

// Code in this file largely adapted from Mercury.  

double improve_mean_to_eccentric_annomaly_guess(const double e, const double M, const double x)
{
  // Based on Mercury
  //  double sx, cx;
  //  sincos(x,&sx,&cx);
      double sx = sin(x);
      double cx = cos(x);
      double es = e*sx;
      double ec = e*cx;
      double f = x-es-M;
      double fp  = 1.-ec;
      double fpp = es;
      double fppp = ec;
      double dx = -f/fp;
      dx = -f/(fp+dx*fpp/2.);
      dx = -f/(fp+dx*fpp/2.+dx*dx*fppp/6.);
      return x+dx;
};

/// calculate mean anomaly to eccentric anomaly
double mean_to_eccentric_annomaly(const double e,  double M)
{
  // Based on Mercury
  const int ORBEL_EHIE_NMAX = 3;

  int nper = static_cast<int>(M/(2.*M_PI));
  M = M - nper*2.*M_PI;
  if(M<0.)  M = M + 2.*M_PI;
  assert(M>=0.);
  if(M>M_PI) 
    { 
      M = 2.*M_PI - M;  
      double x = pow(6.*M,1./3.);
      for(int i=1;i<=ORBEL_EHIE_NMAX;++i)
	x = improve_mean_to_eccentric_annomaly_guess(e,M,x);
      x = 2.*M_PI-x;  
      return x;
    }
  else
    {
      double x = pow(6.*M,1./3.);
      for(int i=1;i<=ORBEL_EHIE_NMAX;++i)
	x = improve_mean_to_eccentric_annomaly_guess(e,M,x);
      return x;
    }
}

/// Calculate the caresian coordinates for the ellipse
void calc_cartesian_for_ellipse(double& x,double& y, double & z, double &vx, double &vy, double &vz, const double a, const double e, const double i, const double O, const double w, const double M, const double GM)
{
  double cape = mean_to_eccentric_annomaly(e,M);
  double scap, ccap;
  sincos(cape,&scap,&ccap);
  double sqe = sqrt(1.-e*e);
  double sqgma = sqrt(GM*a);
  double xfac1 = a*(ccap-e);
  double xfac2 = a*sqe*scap;
  double ri = 1./(a*(1.-e*ccap));
  double vfac1 = -ri*sqgma * scap;
  double vfac2 =  ri*sqgma*sqe*ccap;

  double sw, cw, so, co, si, ci;
  sincos(w,&sw,&cw);
  sincos(O,&so,&co);
  sincos(i,&si,&ci);
  double d1[] = { cw*co-sw*so*ci, cw*so+sw*co*ci, sw*si};
  double d2[] = {-sw*co-cw*so*ci,-sw*so+cw*co*ci, cw*si};
  x  = d1[0]*xfac1+d2[0]*xfac2;
  y  = d1[1]*xfac1+d2[1]*xfac2;
  z  = d1[2]*xfac1+d2[2]*xfac2;
  vx = d1[0]*vfac1+d2[0]*vfac2;
  vy = d1[1]*vfac1+d2[1]*vfac2;
  vz = d1[2]*vfac1+d2[2]*vfac2;
}
void calc_keplerian_for_cartesian( double& a,  double& e,  double& i,  double& O,  double& w,  double& M, const double x,const double y, const double z, const double vx, const double vy, const double vz, const double GM)
{
  const double TINY = 1.e-8;

  double h[] = {y*vz-z*vy, z*vx-x*vz, x*vy-y*vx};
  double h2 = h[0]*h[0]+h[1]*h[1]+h[2]*h[2];
  double hh = sqrt(h2);
  i = acos(h[2]/hh);
  double fac = sqrt(h[0]*h[0]+h[1]*h[1])/hh;
  double u;
  if(fac<TINY)
    {
      O = 0.;
      u = atan2(y,x);
      if(fabs(i-M_PI)<10.*TINY) u = -u;
    }
  else
    {
      O = atan2(h[0],-h[1]);
      u = atan2(z/sin(i),x*cos(O)+y*sin(O));
    }
  if(O<0.) O += 2.*M_PI;
  if(u<0.) u += 2.*M_PI;
  double r = sqrt(x*x+y*y+z*z);
  double energy = (vx*vx+vy*vy+vz*vz)*0.5-GM/r;
  
  if(fabs(energy*r/GM)<sqrt(TINY))
    { // Parabola
      a = 0.5*h2/GM;
      e = 1.;
      double ww = acos(2.*a/r-1.);
      if(vx*x+vy*y+vz*z<0.) w = 2.*M_PI-w;
      double tmpf = tan(0.5*w);
      M = tmpf*(1.+tmpf*tmpf/3.);
      w = u-ww;
      if(w<0.) w+= 2.*M_PI;
      w -= round(w/(2.*M_PI))*2.*M_PI;
    }
  else if (energy<0)
    { // Elipse
      a = -0.5*GM/energy;
      fac = 1.-h2/(GM*a);
      double ww, cape;
      if(fac>TINY)
	{
	  e = sqrt(fac);
	  double face = (a-r)/(a*e);
	  if(face>1.) cape = 0.;
	  else if (face>-1.) cape = acos(face);
	  else cape = M_PI;
	  
	  if(vx*x+vy*y+vz*z<0.) cape = 2.*M_PI-cape;
	  double cw = (cos(cape)-e)/(1.-e*cos(cape));
	  double sw = sqrt(1.-e*e)*sin(cape)/(1.-e*cos(cape));
	  ww = atan2(sw,cw);
	  if(ww<0.) ww += 2.*M_PI;
	}
      else
	{
	  e = 0.;
	  ww = u;
	  cape = u;
	}
      M = cape - e*sin(cape);
      w = u - ww;
      if(w<0.) w += 2.*M_PI;
      w -= round(w/(2.*M_PI))*2.*M_PI;
    }
  else if (energy>0)
    { // Hyperbola
      a = 0.5*GM/energy;
      fac = h2/(GM*a);
      double ww, capf;
      if(fac>TINY)
	{
	  e = sqrt(1.+fac);
	  double tmpf = (a+r)/(a*e);
	  capf = log(tmpf+sqrt(tmpf*tmpf-1.));
	  if(vx*x+vy*y+vz*z<0.) capf = -capf;
	  double cw = (e-cosh(capf))/(e*cosh(capf)-1.);
	  double sw = sqrt(e*e-1.)*sinh(capf)/(e*cosh(capf)-1.);
	  ww = atan2(sw,cw);
	  if(ww<0.) ww += 2.*M_PI;	  
	}
      else
	{
	  e = 1.;
	  double tmpf = 0.5*h2/GM;
	  ww = acos(2.*tmpf/r-1.);
	  if(vx*x+vy*y+vz*z<0.) ww = 2.*M_PI-ww;
	  tmpf = (a+r)/(a*e);
	  capf = log(tmpf+sqrt(tmpf*tmpf-1.));
	}
      M = e*sinh(capf)-capf;
      w = u - ww;
      if(w<0.) w+=2.*M_PI;
      w -= round(w/(2.*M_PI))*2.*M_PI;
    }
};


