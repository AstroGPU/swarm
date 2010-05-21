/*************************************************************************
 * Copyright (C) 2010 by Young In Yeo and the Swarm-NG Development Team  *
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

/*! \file verlet_cpu.cpp
 *  \brief contains member functions of cpu_verlet_integrator and factory to create_cpu_verlet
*/

#include "swarm.h"
#include "verlet_cpu.h"
#include <cassert>
#include <cmath>

namespace swarm {

// factory
extern "C" integrator *create_cpu_verlet(const config &cfg)
{
	return new cpu_verlet_integrator(cfg);
}

/*!
 * \brief Constructor for verlet cpu integrator
 *
 * @param[in] cfg configuration class needs a timestep.
 */
cpu_verlet_integrator::cpu_verlet_integrator(const config &cfg) : 
  m_nbod(0), m_nsys(0), m_is_old_good(0)
{
	if(!cfg.count("time step factor")) ERROR("Integrator cpu_verlet needs a timestep ('time step factor' keyword in the config file).");
	h = atof(cfg.at("time step factor").c_str());
}

/*!
 * \brief function to (re)allocate the integration state
 * 
 * @param[in] ens cpu_ensemble
 */
void cpu_verlet_integrator::alloc_state(cpu_ensemble &ens)
{
	// allocate memory where we'll keep the integration state
	m_nsys = ens.nsys();
	m_nbod = ens.nbod();
	const int len = 3*m_nsys*m_nbod;
	m_xyz_old.resize(len);
	m_vxyz_old.resize(len);
	m_acc.resize(len); m_acc_old.resize(len);
	m_is_old_good = 0;
}


/*!
 * \brief Makes a copy for old data 
 *
 * @param[in,out] ens cpu_ensemble 
 * @param[in] sys system 
 */
void cpu_verlet_integrator::CopyToOld(cpu_ensemble &ens, const unsigned int sys)
{
  for(unsigned int i=0;i<ens.nbod();++i)
    {
      x_old(sys,i)  =  ens.x(sys,i);
      y_old(sys,i)  =  ens.y(sys,i);
      z_old(sys,i)  =  ens.z(sys,i);
      vx_old(sys,i) = ens.vx(sys,i);
      vy_old(sys,i) = ens.vy(sys,i);
      vz_old(sys,i) = ens.vz(sys,i);
      ax_old(sys,i) = ax(sys,i);
      ay_old(sys,i) = ay(sys,i);
      az_old(sys,i) = az(sys,i);
    }
};


/*!
 * \brief Drift step  
 *
 * @param[in,out] ens cpu_ensemble 
 * @param[in] sys system 
 * @param[in] dxdt drift (velocity)  
 * @param[in] dt time step  
 */
void cpu_verlet_integrator::Step_Pos(cpu_ensemble &ens, const unsigned int sys, const std::vector<ThreeVector<real_time> >& dxdt, const real_time dt)
{
	for(unsigned int i=0;i<ens.nbod();++i)
	{
		ens.x(sys,i) += dxdt[i].X() * dt;
		ens.y(sys,i) += dxdt[i].Y() * dt;
		ens.z(sys,i) += dxdt[i].Z() * dt;
	}
};

/*!
 * \brief Kick step  
 *
 * @param[in,out] ens cpu_ensemble 
 * @param[in] sys system 
 * @param[in] dxdt kick  (acceleration)
 * @param[in] dt time step  
 */
void cpu_verlet_integrator::Step_Vel(cpu_ensemble &ens, const unsigned int sys, const std::vector<ThreeVector<real_time> >& dxdt, const real_time dt)
{
	for(unsigned int i=0;i<ens.nbod();++i)
	{
		ens.vx(sys,i) += dxdt[i].X() * dt;
		ens.vy(sys,i) += dxdt[i].Y() * dt;
		ens.vz(sys,i) += dxdt[i].Z() * dt;
	}
};

/*!
 * \brief Updating the time step  
 *
 * @param[in,out] ens cpu_ensemble 
 * @param[in] sys system 
 * @return new time step   
 */
real_time cpu_verlet_integrator::CalcTimeScaleFactor(cpu_ensemble &ens, const unsigned int sys) const
{
	real_time InvPSqEff = 0.;

	for(unsigned int i=0;i<ens.nbod();++i)
	{
		ThreeVector<real_time> xi( ens.x(sys,i), ens.y(sys,i), ens.z(sys,i));
		double msum=ens.mass(sys,i);
		for(unsigned int j=i+1;j<ens.nbod();++j)
		{
			msum += ens.mass(sys,j);
		        ThreeVector<real_time> dist = xi - ThreeVector<real_time>( ens.x(sys,j), ens.y(sys,j), ens.z(sys,j));
			real_time distmagsq = dist.MagnitudeSquared();
			real_time rinv = 1./sqrt(distmagsq);
			rinv *=msum;
			InvPSqEff += rinv/distmagsq;
		}
	}
	return 1./(1.+sqrt(InvPSqEff));  
	// Technically, missing a factor of 2.*M_PI, but this choice is arbitrary 
};

/*!
 * \brief Calculation for the drift  
 *
 * @param[in] ens cpu_ensemble 
 * @param[in] sys system 
 * @return drift (velocity)  
 */
std::vector<ThreeVector<real_time> > cpu_verlet_integrator::CalcDerivForDrift(cpu_ensemble &ens, const unsigned int sys) const
{
	
	std::vector<ThreeVector<real_time> > vi(ens.nbod(),ThreeVector<real_time>(0.,0.,0.));
	
	for(unsigned int i=0;i<ens.nbod();++i)
	{
		vi[i].X()=ens.vx(sys,i);
		vi[i].Y()=ens.vy(sys,i);
		vi[i].Z()=ens.vz(sys,i);
	}
	return vi;
};

/*!
 * \brief Calculation for the kick  
 *
 * @param[in] ens cpu_ensemble 
 * @param[in] sys system 
 * @return kick (acceleration) 
 */
std::vector<ThreeVector<real_time> > cpu_verlet_integrator::CalcDerivForKick(cpu_ensemble &ens, const unsigned int sys) const
{
	std::vector<ThreeVector<real_time> > acc(ens.nbod(),ThreeVector<real_time>(0.,0.,0.));

	for(unsigned int i=1;i<ens.nbod();++i)
	{
		ThreeVector<real_time> xi( ens.x(sys,i), ens.y(sys,i), ens.z(sys,i));
		for(unsigned int j=i+1;j<ens.nbod();++j)
		{
		        ThreeVector<real_time> dist = xi - ThreeVector<real_time>( ens.x(sys,j), ens.y(sys,j), ens.z(sys,j));
			real_time distmagsq = dist.MagnitudeSquared();
			real_time distmagcube = distmagsq*sqrt(distmagsq);
			acc[i].X() -= ens.mass(sys,j)*dist.X()/distmagcube;
			acc[i].Y() -= ens.mass(sys,j)*dist.Y()/distmagcube;
			acc[i].Z() -= ens.mass(sys,j)*dist.Z()/distmagcube;
			acc[j].X() += ens.mass(sys,i)*dist.X()/distmagcube;
			acc[j].Y() += ens.mass(sys,i)*dist.Y()/distmagcube;
			acc[j].Z() += ens.mass(sys,i)*dist.Z()/distmagcube;
		}
	}
	// Add in star-planet accel last to reduce roundoff error
	{
		unsigned int i=0;
		ThreeVector<real_time> xi( ens.x(sys,i), ens.y(sys,i), ens.z(sys,i));
		for(unsigned int j=1;j<ens.nbod();++j)
		{
		        ThreeVector<real_time> dist = xi - ThreeVector<real_time>( ens.x(sys,j), ens.y(sys,j), ens.z(sys,j));
			real_time distmagsq = dist.MagnitudeSquared();
			real_time distmagcube = distmagsq*sqrt(distmagsq);
			acc[i].X() -= ens.mass(sys,j)*dist.X()/distmagcube;
			acc[i].Y() -= ens.mass(sys,j)*dist.Y()/distmagcube;
			acc[i].Z() -= ens.mass(sys,j)*dist.Z()/distmagcube;
			acc[j].X() += ens.mass(sys,i)*dist.X()/distmagcube;
			acc[j].Y() += ens.mass(sys,i)*dist.Y()/distmagcube;
			acc[j].Z() += ens.mass(sys,i)*dist.Z()/distmagcube;
		}
	}
	return acc;
};

/*!
 * \brief Verlet CPU integrate function 
 *
 * Verlet GPU implementation is based on this
 * @see swarm::prop_verlet 
 * @param[in,out] ens cpu_ensemble 
 * @param[in] dT destination time 
 */
void cpu_verlet_integrator::integrate(cpu_ensemble &ens, real_time dT)
{
	// Allocate integration state (if not already there)
	if(ens.last_integrator() != this) 
	  {
	    alloc_state(ens);
	    ens.set_last_integrator(this);
	    if(dT == 0.) { return; }
	  }


	for ( unsigned int sys=0;sys<ens.nsys();++sys )
	{

	  if(!ens.is_active(sys)) continue;

	  real_time Tend = ens.time( sys ) + dT;

	  std::vector<ThreeVector<real_time> > dposdt(ens.nbod()), dveldt(ens.nbod());
	  // propagate the system until we match or exceed Tend
	  while ( ens.time( sys ) < Tend )
	    {
	  double T = ens.time(sys);
	  //double h = T + this->h <= Tend ? this->h : Tend - T;
	  double h_half= h*0.5;

	  //real_time HalfDeltaTLast = h*CalcTimeScaleFactor(ens,sys);
	  real_time HalfDeltaTLast = h_half;
		    h_half = HalfDeltaTLast;
		    dposdt = CalcDerivForDrift(ens,sys);
		    Step_Pos(ens,sys,dposdt,h_half);

		    dveldt = CalcDerivForKick(ens,sys);	
		    Step_Vel(ens,sys,dveldt,h_half);

	      ens.time( sys ) += h_half;
	//	    Time += h_half;

		    h_half =1./(2./(CalcTimeScaleFactor(ens,sys)*h)-1./h_half);
		    //h_half =(CalcTimeScaleFactor(ens,sys));
		    Step_Vel(ens,sys,dveldt,h_half);	

		    dposdt = CalcDerivForDrift(ens,sys);
		    Step_Pos(ens,sys,dposdt,h_half);

		    //Time += h_half;
		    HalfDeltaTLast = h_half;
		    //++mnSteps;

//	      Evolve ( ens,sys );
	      ens.time( sys ) += h_half;
	      //ens.time( sys ) += h;
	      ens.nstep(sys)++;
	    }
	} // end loop over systems
}

} // end namespace swarm
