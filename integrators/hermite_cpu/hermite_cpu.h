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

/*! \file hermite_cpu.h
 *  \brief declares cpu_hermite_integrator
 *
 *  Note that while this clas derivers from integrator, it does not use cpu_generic_integrator
*/


#ifndef integ_hermite_cpu_h__
#define integ_hermite_cpu_h__

#include "swarm.h"
#include <valarray>

namespace swarm {

/*!
 * \brief cpu_hermite_integrator class
 */
class cpu_hermite_integrator : public integrator
{
public:
	typedef double real;
	typedef real   real_time;
	typedef real   real_mass;
	typedef real   real_pos;
	typedef real   real_vel;
	typedef real   real_acc;
	typedef real   real_jerk;

protected:
	real_time m_h;

	// Integration state
	std::valarray<real_pos>		m_xyz_old;
	std::valarray<real_vel>		m_vxyz_old;
	std::valarray<real_acc>		m_acc, m_acc_old;
	std::valarray<real_jerk>	m_jerk, m_jerk_old;
	int				m_nsys, m_nbod;
	int                             m_is_old_good;

	/*!
  	 * \brief  function to (re)allocate the integration state
	 * 
	 * @param[in] ens cpu_ensemble
 	 */
	void alloc_state(cpu_ensemble &ens);

public:
	cpu_hermite_integrator(const config &cfg);

public:
	virtual void integrate(cpu_ensemble &ens, real_time T);

	// Is it dangerous to provide these as public?
	// If not, people could use them to interpolate to some specific time
	void set_timestep(const real_time hnew) { m_h = hnew; };

	int is_old_good() const { return m_is_old_good; };
	
protected:
	void predict(cpu_ensemble &ens, const unsigned int sys, real_time hh);
	void Correct(cpu_ensemble &ens, const unsigned int sys, real_time hh);
	void CorrectAlpha7by6(cpu_ensemble &ens, const unsigned int sys, real_time hh);

	void UpdateAccJerk(cpu_ensemble &ens, const unsigned int sys);
	
	void Evolve(cpu_ensemble &ens, const unsigned int sys, real_time hh);
	void EvolvePEC1(cpu_ensemble &ens, const unsigned int sys, real_time hh);
	void EvolvePEC2(cpu_ensemble &ens, const unsigned int sys, real_time hh);
 
	void CopyToOld(cpu_ensemble &ens, const unsigned int sys);

	public:
	// non-const versions
	/// returns reference to x component of acceleration for system sys and body bod
	double&  ax(int sys, int bod) { return m_acc[bod*m_nsys + sys]; }
	double&  ay(int sys, int bod) { return m_acc[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double&  az(int sys, int bod) { return m_acc[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	double& jx(int sys, int bod) { return m_jerk[bod*m_nsys + sys]; }
	double& jy(int sys, int bod) { return m_jerk[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double& jz(int sys, int bod) { return m_jerk[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	double&  x_old(int sys, int bod) { return m_xyz_old[bod*m_nsys + sys]; }
	double&  y_old(int sys, int bod) { return m_xyz_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double&  z_old(int sys, int bod) { return m_xyz_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	double& vx_old(int sys, int bod) { return m_vxyz_old[bod*m_nsys + sys]; }
	double& vy_old(int sys, int bod) { return m_vxyz_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double& vz_old(int sys, int bod) { return m_vxyz_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	double&  ax_old(int sys, int bod) { return m_acc_old[bod*m_nsys + sys]; }
	double&  ay_old(int sys, int bod) { return m_acc_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double&  az_old(int sys, int bod) { return m_acc_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	double& jx_old(int sys, int bod) { return m_jerk_old[bod*m_nsys + sys]; }
	double& jy_old(int sys, int bod) { return m_jerk_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double& jz_old(int sys, int bod) { return m_jerk_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	// const versions
	double  ax(int sys, int bod) const { return m_acc[bod*m_nsys + sys]; }
	double  ay(int sys, int bod) const { return m_acc[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double  az(int sys, int bod) const { return m_acc[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	double jx(int sys, int bod) const { return m_jerk[bod*m_nsys + sys]; }
	double jy(int sys, int bod) const { return m_jerk[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double jz(int sys, int bod) const { return m_jerk[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	double  x_old(int sys, int bod) const { return m_xyz_old[bod*m_nsys + sys]; }
	double  y_old(int sys, int bod) const { return m_xyz_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double  z_old(int sys, int bod) const { return m_xyz_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	double vx_old(int sys, int bod) const { return m_vxyz_old[bod*m_nsys + sys]; }
	double vy_old(int sys, int bod) const { return m_vxyz_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double vz_old(int sys, int bod) const { return m_vxyz_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	double  ax_old(int sys, int bod) const { return m_acc_old[bod*m_nsys + sys]; }
	double  ay_old(int sys, int bod) const { return m_acc_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double  az_old(int sys, int bod) const { return m_acc_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	double jx_old(int sys, int bod) const { return m_jerk_old[bod*m_nsys + sys]; }
	double jy_old(int sys, int bod) const { return m_jerk_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double jz_old(int sys, int bod) const { return m_jerk_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }
};

} // end namespace swarm
#endif
