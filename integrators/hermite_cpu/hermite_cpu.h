#ifndef integ_hermite_cpu_h__
#define integ_hermite_cpu_h__

#include "ThreeVector.hpp"
#include "swarm.h"
#include <valarray>

#if 0  
// Eventually, we'll combine my hermite_cpu with the class that Young
// In writes for his gpu-based hermite integrator.  I hope  my code
// can be integrated by searching for hermite_cpu and replacing it
// with hermite.

class gpu_hermite_cpu_integrator : public integrator
{
 protected:
	float h;
	
	dim3 gridDim;
	int threadsPerBlock;

public:
	gpu_hermite_cpu_integrator(const config &cfg);

public:
	// Young In will put his code into a similar place 
	// but initially in a class with a different name
	// eventually we'll merge these two, so his code would go here
	void integrate(gpu_ensemble &ens, float T)  { abort(); }

	// No support for CPU execution. Note: we could make this function
	// transparently copy the CPU ensemble to the GPU, and back once
	// the integration is done
	void integrate(cpu_ensemble &ens, float T) { abort(); }
};
#endif



class cpu_hermite_cpu_integrator : public integrator
{
public:
	typedef double real;
	typedef float  real_time;
	typedef real   real_mass;
	typedef real   real_pos;
	typedef real   real_vel;
	typedef real   real_acc;
	typedef real   real_jerk;

protected:
	float h;

	// Integration state
	std::valarray<real_pos>		m_xyz_old;
	std::valarray<real_vel>		m_vxyz_old;
	std::valarray<real_acc>		m_acc, m_acc_old;
	std::valarray<real_jerk>	m_jerk, m_jerk_old;
	int				m_nsys, m_nbod;
	void alloc_state(cpu_ensemble &ens);

public:

	cpu_hermite_cpu_integrator(const config &cfg);
//	~cpu_hermite_cpu_integrator();

public:
	// No support for GPU execution. Note: we could make this function
	// transparently copy the CPU ensemble to the GPU, and back once
	// the integration is done
	virtual void integrate(gpu_ensemble &ens, real_time T) { abort(); }
	virtual void integrate(cpu_ensemble &ens, real_time T);

	void predict(cpu_ensemble &ens, const unsigned int sys);
	

// 	virtual void allocate_workspace_for_ensemble(const ensemble &ens)
// 	{
// 	  if(m_nsys != ens.nsys() || m_nbod != ens.nbod())
// 	    {
// 
// 	      m_xyz_old = (double*)realloc(m_xyz_old, 3*m_nsys*m_nbod*sizeof(*m_xyz_old));
// 	      m_vxyz_old = (double*)realloc(m_vxyz_old, 3*m_nsys*m_nbod*sizeof(*m_vxyz_old));
// 	      m_acc = (double*)realloc(m_acc, 3*m_nsys*m_nbod*sizeof(*m_acc));
// 	      m_jerk = (double*)realloc(m_jerk, 3*m_nsys*m_nbod*sizeof(*m_jerk));
// 	      m_acc_old = (double*)realloc(m_acc_old, 3*m_nsys*m_nbod*sizeof(*m_acc_old));
// 	      m_jerk_old = (double*)realloc(m_jerk_old, 3*m_nsys*m_nbod*sizeof(*m_jerk_old));
// 	      
// 	      if(m_xyz_old==NULL) ERROR("CPU memory allocation failed!\n");
// 	      if(m_vxyz_old==NULL) ERROR("CPU memory allocation failed!\n");
// 	      if(m_acc==NULL) ERROR("CPU memory allocation failed!\n");
// 	      if(m_jerk==NULL) ERROR("CPU memory allocation failed!\n");
// 	      if(m_acc_old==NULL) ERROR("CPU memory allocation failed!\n");
// 	      if(m_jerk_old==NULL) ERROR("CPU memory allocation failed!\n");
// 
// 	      m_nsys = ens.nsys();
// 	      m_nbod = ens.nbod();
// 	    }
// 	}

/*	void deallocate_workspace()
	{

	  fprintf(stderr, "Warning: Deallocating workspace for cpu_hermite\n");
	  ::free(m_xyz_old); m_xyz_old = NULL;
	  ::free(m_vxyz_old);m_vxyz_old = NULL;
	  ::free(m_acc); m_acc = NULL;
	  ::free(m_jerk); m_jerk = NULL;
	  ::free(m_acc_old); m_acc_old = NULL;
	  ::free(m_jerk_old); m_jerk_old = NULL;
	}
*/
 protected:
	void Correct(cpu_ensemble &ens, const unsigned int sys);
	void CorrectAlpha7by6(cpu_ensemble &ens, const unsigned int sys);

	void UpdateAccJerk(cpu_ensemble &ens, const unsigned int sys);
	
	void Evolve(cpu_ensemble &ens, const unsigned int sys);
	void EvolvePEC1(cpu_ensemble &ens, const unsigned int sys);
	void EvolvePEC2(cpu_ensemble &ens, const unsigned int sys);
  
	void CopyToOld(cpu_ensemble &ens, const unsigned int sys);

	public:
	// non-const versions
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

#endif
