#ifndef integ_verlet_cpu_h__
#define integ_verlet_cpu_h__

#include "swarm.h"
#include "ThreeVector.hpp"
#include <valarray>
#include <vector>

namespace swarm {

/**
 * cpu_verlet_integrator class
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

	/**
	 * function to (re)allocate the integration state
	 * 
	 * @param[in] ens cpu_ensemble
	 */
	void alloc_state(cpu_ensemble &ens);

public:
	/**
	 * Constructor for verlet cpu integrator
	 *
	 * @param[in] cfg configuration file needs a timestep.
	 */
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
	double&  ax(int sys, int bod) { return m_acc[bod*m_nsys + sys]; }
	double&  ay(int sys, int bod) { return m_acc[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double&  az(int sys, int bod) { return m_acc[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }


	double&  x_old(int sys, int bod) { return m_xyz_old[bod*m_nsys + sys]; }
	double&  y_old(int sys, int bod) { return m_xyz_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double&  z_old(int sys, int bod) { return m_xyz_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	double& vx_old(int sys, int bod) { return m_vxyz_old[bod*m_nsys + sys]; }
	double& vy_old(int sys, int bod) { return m_vxyz_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double& vz_old(int sys, int bod) { return m_vxyz_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	double&  ax_old(int sys, int bod) { return m_acc_old[bod*m_nsys + sys]; }
	double&  ay_old(int sys, int bod) { return m_acc_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double&  az_old(int sys, int bod) { return m_acc_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }


	// const versions
	double  ax(int sys, int bod) const { return m_acc[bod*m_nsys + sys]; }
	double  ay(int sys, int bod) const { return m_acc[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double  az(int sys, int bod) const { return m_acc[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }


	double  x_old(int sys, int bod) const { return m_xyz_old[bod*m_nsys + sys]; }
	double  y_old(int sys, int bod) const { return m_xyz_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double  z_old(int sys, int bod) const { return m_xyz_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	double vx_old(int sys, int bod) const { return m_vxyz_old[bod*m_nsys + sys]; }
	double vy_old(int sys, int bod) const { return m_vxyz_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double vz_old(int sys, int bod) const { return m_vxyz_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	double  ax_old(int sys, int bod) const { return m_acc_old[bod*m_nsys + sys]; }
	double  ay_old(int sys, int bod) const { return m_acc_old[m_nbod*m_nsys + bod*m_nsys + sys]; }
	double  az_old(int sys, int bod) const { return m_acc_old[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

	void Step_Pos(cpu_ensemble &ens, const unsigned int sys, const std::vector<ThreeVector<real> >& dxdt, const real dt);
	void Step_Vel(cpu_ensemble &ens, const unsigned int sys, const std::vector<ThreeVector<real> >& dxdt, const real dt);
	real CalcTimeScaleFactor(cpu_ensemble &ens, const unsigned int sys) const;
	std::vector<ThreeVector<real> > CalcDerivForDrift(cpu_ensemble &ens, const unsigned int sys) const;
	std::vector<ThreeVector<real> > CalcDerivForKick(cpu_ensemble &ens, const unsigned int sys) const;


};

} // end namespace swarm
#endif
