#ifndef integ_hermite_adap_h__
#define integ_hermite_adap_h__

#include "swarm.h"

#define LARGE_NUMBER 1e16
#define SMALL_NUMBER 1e-16
#define STEP_FACTOR 0.01
#define NBODIES 3

typedef double real;
typedef real   real_time;
typedef float   real_mass;
typedef real   real_pos;
typedef real   real_vel;
typedef real   real_acc;
typedef real   real_jerk;

struct gpu_hermite_adap_integrator_data
{
       // Integration state
       real_pos        *m_xyz_old;
       real_vel        *m_vxyz_old;
       real_acc        *m_acc, *m_acc_old;
       real_jerk       *m_jerk, *m_jerk_old;

       // Other parameters the integrator will need
       // For the adaptive time step, h is now the minimum
       float h;
};


//class gpu_hermite_integrator : public integrator
class gpu_hermite_adap_integrator : public integrator, public gpu_hermite_adap_integrator_data
{
protected:
	float h;
	int prec;
	
	dim3 gridDim;
	int threadsPerBlock;

public:
	gpu_hermite_adap_integrator(const config &cfg);

public:
	void integrate(gpu_ensemble &ens, double T);
};

#endif
