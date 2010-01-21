#ifndef integ_hermite_h__
#define integ_hermite_h__

#include "swarm.h"

typedef double real;
typedef real   real_time;
typedef real   real_mass;
typedef real   real_pos;
typedef real   real_vel;
typedef real   real_acc;
typedef real   real_jerk;

struct gpu_hermite_integrator_data
{
       // Integration state
       real_pos        *m_xyz_old;
       real_vel        *m_vxyz_old;
       real_acc        *m_acc, *m_acc_old;
       real_jerk       *m_jerk, *m_jerk_old;

       // Other parameters the integrator will need
       float h;
};


//class gpu_hermite_integrator : public integrator
class gpu_hermite_integrator : public integrator, public gpu_hermite_integrator_data
{
protected:
	float h;
	
	dim3 gridDim;
	int threadsPerBlock;

public:
	gpu_hermite_integrator(const config &cfg);

public:
	void integrate(gpu_ensemble &ens, double T);
};

#endif
