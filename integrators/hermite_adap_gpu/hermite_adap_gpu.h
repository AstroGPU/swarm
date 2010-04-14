#ifndef integ_hermite_adap_h__
#define integ_hermite_adap_h__

#include "swarm.h"
#include "user.h"

#define LARGE_NUMBER 1e16
#define SMALL_NUMBER 1e-16

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
       // h is the minimum time step allowed
       // stepfac is the fraction of the implicit step you wish to take, e.g., dt = acc.mag/jerk.mag * stepfac + h
       float h, stepfac;
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
