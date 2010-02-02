#ifndef integ_hermite_h__
#define integ_hermite_h__

#include "swarm.h"

namespace swarm {
namespace hermite_gpu {

template<typename real_hi, typename real_lo>
struct gpu_hermite_integrator_data
{
  // Typedef's (not clear if Jianwei is using these throughout
  typedef real_hi real;
  typedef real_hi real_time;
  typedef float   real_mass;
  typedef real_hi real_pos;
  typedef real_hi real_vel;
  typedef real_lo real_acc;
  typedef real_lo real_jerk;

  // Integration state
  real_pos        *m_xyz_old;
  real_vel        *m_vxyz_old;
  real_acc        *m_acc, *m_acc_old;
  real_jerk       *m_jerk, *m_jerk_old;
  
  // Other parameters the integrator will need
  float h;
};


//class gpu_hermite_integrator : public integrator
template< typename real_hi, typename real_lo>
class gpu_hermite_integrator : public integrator, public gpu_hermite_integrator_data<real_hi,real_lo>
{
protected:
	float h;
	int prec;
	
	dim3 gridDim;
	int threadsPerBlock;

public:
	gpu_hermite_integrator(const config &cfg);

public:
	void integrate(gpu_ensemble &ens, double T);

};




} // end namespace hermite_gpu
} // end namespace swarm

#endif
