#include "swarm.h"
#include "hermite_cpu.h"
#include "ThreeVector.hpp"
#include <cassert>
#include <cmath>

// We don't need to use this for CPU integrators
// ensemble cpu_hermite_ens;

// factory
extern "C" integrator *create_cpu_hermite(const config &cfg)
{
	return new cpu_hermite_integrator(cfg);
}

// Constructor
cpu_hermite_integrator::cpu_hermite_integrator(const config &cfg) : 
  m_nbod(0), m_nsys(0), m_is_old_good(0)
{
	if(!cfg.count("h")) ERROR("Integrator cpu_hermite needs a timestep ('h' keyword in the config file).");
	h = atof(cfg.at("h").c_str());
}

void cpu_hermite_integrator::alloc_state(cpu_ensemble &ens)
{
	// allocate memory where we'll keep the integration state
	m_nsys = ens.nsys();
	m_nbod = ens.nbod();
	const int len = 3*m_nsys*m_nbod;
	m_xyz_old.resize(len);
	m_vxyz_old.resize(len);
	m_acc.resize(len); m_acc_old.resize(len);
	m_jerk.resize(len); m_jerk_old.resize(len);
	m_is_old_good = 0;
}

// Actual calculations
void cpu_hermite_integrator::predict(cpu_ensemble &ens, const unsigned int sys)
{
  for(unsigned int i=0;i<ens.nbod();++i)
    {
      const real_time hby2 = h/2., hby3 = h/3.;
      ens.x(sys,i)  += h*(ens.vx(sys,i)+hby2*(ax(sys,i)+hby3*jx(sys,i)));
      ens.y(sys,i)  += h*(ens.vy(sys,i)+hby2*(ay(sys,i)+hby3*jy(sys,i)));
      ens.z(sys,i)  += h*(ens.vz(sys,i)+hby2*(az(sys,i)+hby3*jz(sys,i)));
      ens.vx(sys,i) += h*(ax(sys,i)+hby2*jx(sys,i));
      ens.vy(sys,i) += h*(ay(sys,i)+hby2*jy(sys,i));
      ens.vz(sys,i) += h*(az(sys,i)+hby2*jz(sys,i));
    }
};

void cpu_hermite_integrator::CopyToOld(cpu_ensemble &ens, const unsigned int sys)
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
      jx_old(sys,i) = jx(sys,i);
      jy_old(sys,i) = jy(sys,i);
      jz_old(sys,i) = jz(sys,i);
    }
};

void cpu_hermite_integrator::CorrectAlpha7by6(cpu_ensemble &ens, const unsigned int sys)
{
  for(unsigned int i=0;i<ens.nbod();++i)
    {
      const real_time hby2 = h/2., hby6 = h/6.;
      const real_time h7by30 = h*7./30, hby7 = h/7.;
      ens.vx(sys,i) = vx_old(sys,i) + hby2*((ax_old(sys,i)+ax(sys,i)
					 + hby6*(jx_old(sys,i)-jx(sys,i))));
      ens.vy(sys,i) = vy_old(sys,i) + hby2*((ay_old(sys,i)+ay(sys,i)
					 + hby6*(jy_old(sys,i)-jy(sys,i))));
      ens.vz(sys,i) = vz_old(sys,i) + hby2*((az_old(sys,i)+az(sys,i)
					 + hby6*(jz_old(sys,i)-jz(sys,i))));
      ens.x(sys,i) = x_old(sys,i) + hby2*((vx_old(sys,i)+ens.vx(sys,i))
				      + h7by30*((ax_old(sys,i)-ax(sys,i))
				      + hby7*(jx_old(sys,i)+jx(sys,i))));
      ens.y(sys,i) = y_old(sys,i) + hby2*((vy_old(sys,i)+ens.vy(sys,i))
				      + h7by30*((ay_old(sys,i)-ay(sys,i))
				      + hby7*(jy_old(sys,i)+jy(sys,i))));
      ens.z(sys,i) = z_old(sys,i) + hby2*((vz_old(sys,i)+ens.vz(sys,i))
				      + h7by30*((az_old(sys,i)-az(sys,i))
				      + hby7*(jz_old(sys,i)+jz(sys,i))));
    }
};


// Kokubo & Makino PASJ 56, 861 explain choice of alpha = 7/6
void cpu_hermite_integrator::Correct(cpu_ensemble &ens, const unsigned int sys)
{
  CorrectAlpha7by6(ens,sys);
}

void cpu_hermite_integrator::UpdateAccJerk(cpu_ensemble &ens, const unsigned int sys)
{
  // Use's double precission
  ThreeVector<real> dx(0.,0.,0.), dv(0.,0.,0.);
  
    { // First calculate acceleration and jerk for the Sun
      unsigned int i = 0;
      ThreeVector<real> xi( ens.x(sys,i), ens.y(sys,i), ens.z(sys,i));
      ThreeVector<real> vi(ens.vx(sys,i),ens.vy(sys,i),ens.vz(sys,i));
      ThreeVector<real> ai(0.,0.,0.), ji(0.,0.,0.);
      for(unsigned int j=1;j<ens.nbod();++j)
	{
	  //	    dx = mPos[j] - mPos[i];
	  dx = ThreeVector<real>( ens.x(sys,j), ens.y(sys,j), ens.z(sys,j)) - xi;
	  //	    dv = mVel[j] - mVel[i];
	  dv = ThreeVector<real>(ens.vx(sys,j),ens.vy(sys,j),ens.vz(sys,j))  - vi;
	  real r2 = dx.MagnitudeSquared();
	  real rv = dot(dx,dv);
	  real rinv = 1./sqrt(r2);
	  rv *= 3./r2;
	  rinv *= ens.mass(sys,j); // mMass[j];
	  real rinv3 = rinv/r2;
	  
	  dx *= rinv3;
	  ai += dx;
	  dv *= rinv3;
	  ji += dv;
	  dx *= rv;
	  ji -= dx;
	}
      ax(sys,i) = ai.X();
      ay(sys,i) = ai.Y();
      az(sys,i) = ai.Z();
      jx(sys,i) = ji.X();
      jy(sys,i) = ji.Y();
      jz(sys,i) = ji.Z();
    }
  
  // Now calculate acceleration and jerk for the planets
  for(unsigned int i=1;i<ens.nbod();++i)
    {
      ThreeVector<real> xi( ens.x(sys,i), ens.y(sys,i), ens.z(sys,i));
      ThreeVector<real> vi(ens.vx(sys,i),ens.vy(sys,i),ens.vz(sys,i));
      ThreeVector<real> ai(0.,0.,0.), ji(0.,0.,0.);
      for(unsigned int j=1;j<ens.nbod();++j)
	{
	  if(j==i) continue; // Ignore body interacting with itself
	  //	    dx = mPos[j] - mPos[i];
	  dx = ThreeVector<real>( ens.x(sys,j), ens.y(sys,j), ens.z(sys,j)) - xi;
	  //	    dv = mVel[j] - mVel[i];
	  dv = ThreeVector<real>(ens.vx(sys,j),ens.vy(sys,j),ens.vz(sys,j)) - vi;
	  real r2 = dx.MagnitudeSquared();
	  real rv = dot(dx,dv);
	  real rinv = 1./sqrt(r2);
	  rv *= 3./r2;
	  rinv *= ens.mass(sys,j); // mMass[j];
	  real rinv3 = rinv/r2;
	  
	  dx *= rinv3;
	  ai += dx;
	  dv *= rinv3;
	  ji += dv;
	  dx *= rv;
	  ji -= dx;
	}
      
      {  // But add sun's contribution last to minimize round-off error
	unsigned int j=0;
	//	    dx = mPos[j] - mPos[i];
	dx = ThreeVector<real>( ens.x(sys,j), ens.y(sys,j), ens.z(sys,j)) - xi;	      
	//	    dv = mVel[j] - mVel[i];
	//	    dv = mVel[j] - vi;	    
	dv = ThreeVector<real>(ens.vx(sys,j),ens.vy(sys,j),ens.vz(sys,j)) - vi;
	real r2 = dx.MagnitudeSquared();
	real rv = dot(dx,dv);
	real rinv = 1./sqrt(r2);
	rv *= 3./r2;
	rinv *= ens.mass(sys,j); // mMass[j];
	real rinv3 = rinv/r2;
	
	dx *= rinv3;
	ai += dx;
	dv *= rinv3;
	ji += dv;
	dx *= rv;
	ji -= dx;
	
      } // end adding Sun's contribution
      ax(sys,i) = ai.X();
      ay(sys,i) = ai.Y();
      az(sys,i) = ai.Z();
      jx(sys,i) = ji.X();
      jy(sys,i) = ji.Y();
      jz(sys,i) = ji.Z();       
    } // end loop over bodies
}


void cpu_hermite_integrator::EvolvePEC2(cpu_ensemble &ens, const unsigned int sys)
{
  CopyToOld(ens,sys);
  predict(ens,sys);
  UpdateAccJerk(ens,sys);
  Correct(ens,sys);
  UpdateAccJerk(ens,sys); 
  Correct(ens,sys);
};

  
// See Kokubo, Yoshinaga & Makino MNRAS 297, 1067 for details of choice
void cpu_hermite_integrator::Evolve(cpu_ensemble &ens, const unsigned int sys)
{
  EvolvePEC2(ens,sys);
  m_is_old_good = 1;
};

void cpu_hermite_integrator::integrate(cpu_ensemble &ens, real_time dT)
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
	  if(!ens.active(sys)) continue;
	  real_time Tend = ens.time( sys ) + dT;
	  UpdateAccJerk ( ens,sys );
	  // propagate the system until we match or exceed Tend
	  while ( ens.time( sys ) < Tend )
	    {
	      Evolve ( ens,sys );
	      ens.time( sys ) += h;
	    }
	} // end loop over systems
}
