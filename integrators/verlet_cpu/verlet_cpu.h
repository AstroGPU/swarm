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

/*! \file verlet_cpu.h
 * \brief declares cpu_verlet_integrator class and inline accessor functions
*/


#ifndef integ_verlet_cpu_h__
#define integ_verlet_cpu_h__

#include "swarm.h"
#include "ThreeVector.hpp"
#include <valarray>
#include <vector>

namespace swarm {

/*!
 * \brief cpu_verlet_integrator class
 */
class cpu_verlet_integrator : public integrator
{
public:
	typedef double real;
	typedef real   real_time;
	typedef real   real_mass;
	typedef real   real_pos;
	typedef real   real_vel;
	typedef real   real_acc;

protected:
	real_time h;

	// Integration state
	std::valarray<real_pos>		m_xyz_old;
	std::valarray<real_vel>		m_vxyz_old;
	std::valarray<real_acc>		m_acc, m_acc_old;
	int				m_nsys, m_nbod;
	int                             m_is_old_good;

	void alloc_state(cpu_ensemble &ens);

public:
	cpu_verlet_integrator(const config &cfg);

public:
	virtual void integrate(cpu_ensemble &ens, real_time T);

	void set_timestep(const real_time hnew) { h = hnew; };

	int is_old_good() const { return m_is_old_good; };
	
protected:
	void predict(cpu_ensemble &ens, const unsigned int sys);

	void CopyToOld(cpu_ensemble &ens, const unsigned int sys);

	public:
	// non-const versions
	/// returns reference to x component of acceleration for system sys and body bod
	double&  ax(int sys, int bod) { return m_acc[bod*m_nsys + sys]; }
	/// returns reference to y component of acceleration for system sys and body bod
	double&  ay(int sys, int bod) { return m_acc[m_nbod*m_nsys + bod*m_nsys + sys]; }
	/// returns reference to z component of acceleration for system sys and body bod
	double&  az(int sys, int bod) { return m_acc[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }


	/// returns reference to x coordinate for system sys and body bod at previous time
	double&  x_old(int sys, int bod) { return m_xyz_old[bod*m_nsys + sys]; }
	/// returns reference to y coordinate for system sys and body bod at previous time
	double&  y_old(int sys, int bod) { return m_xyz_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	/// returns reference to z coordinate for system sys and body bod at previous time
	double&  z_old(int sys, int bod) { return m_xyz_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	/// returns reference to x component of velocity for system sys and body bod at previous time
	double& vx_old(int sys, int bod) { return m_vxyz_old[bod*m_nsys + sys]; }
	/// returns reference to y component of velocity for system sys and body bod at previous time
	double& vy_old(int sys, int bod) { return m_vxyz_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	/// returns reference to z component of velocity for system sys and body bod at previous time
	double& vz_old(int sys, int bod) { return m_vxyz_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	/// returns reference to x component of acceleration for system sys and body bod at previous time
	double&  ax_old(int sys, int bod) { return m_acc_old[bod*m_nsys + sys]; }
	/// returns reference to y component of acceleration for system sys and body bod at previous time
	double&  ay_old(int sys, int bod) { return m_acc_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	/// returns reference to z component of acceleration for system sys and body bod at previous time
	double&  az_old(int sys, int bod) { return m_acc_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }


	// const versions
	/// returns value of x component of acceleration for system sys and body bod
	double  ax(int sys, int bod) const { return m_acc[bod*m_nsys + sys]; }
	/// returns value of y component of acceleration for system sys and body bod
	double  ay(int sys, int bod) const { return m_acc[m_nbod*m_nsys + bod*m_nsys + sys]; }
	/// returns value of z component of acceleration for system sys and body bod
	double  az(int sys, int bod) const { return m_acc[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }


	/// returns value of x coordinate for system sys and body bod at previous time
	double  x_old(int sys, int bod) const { return m_xyz_old[bod*m_nsys + sys]; }
	/// returns value of y coordinate for system sys and body bod at previous time
	double  y_old(int sys, int bod) const { return m_xyz_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	/// returns value of z coordinate for system sys and body bod at previous time
	double  z_old(int sys, int bod) const { return m_xyz_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	/// returns value of x component of velocity for system sys and body bod at previous time
	double vx_old(int sys, int bod) const { return m_vxyz_old[bod*m_nsys + sys]; }
	/// returns value of y component of velocity for system sys and body bod at previous time
	double vy_old(int sys, int bod) const { return m_vxyz_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	/// returns value of z component of velocity for system sys and body bod at previous time
	double vz_old(int sys, int bod) const { return m_vxyz_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	/// returns value of x component of acceleration for system sys and body bod at previous time
	double  ax_old(int sys, int bod) const { return m_acc_old[bod*m_nsys + sys]; }
	/// returns value of y component of acceleration for system sys and body bod at previous time
	double  ay_old(int sys, int bod) const { return m_acc_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	/// returns value of z component of acceleration for system sys and body bod at previous time
	double  az_old(int sys, int bod) const { return m_acc_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	void Step_Pos(cpu_ensemble &ens, const unsigned int sys, const std::vector<ThreeVector<real> >& dxdt, const real dt);
	void Step_Vel(cpu_ensemble &ens, const unsigned int sys, const std::vector<ThreeVector<real> >& dxdt, const real dt);
	real CalcTimeScaleFactor(cpu_ensemble &ens, const unsigned int sys) const;
	std::vector<ThreeVector<real> > CalcDerivForDrift(cpu_ensemble &ens, const unsigned int sys) const;
	std::vector<ThreeVector<real> > CalcDerivForKick(cpu_ensemble &ens, const unsigned int sys) const;


};

} // end namespace swarm
#endif
