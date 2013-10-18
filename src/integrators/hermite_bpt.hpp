/*************************************************************************
 * Copyright (C) 2011 by Saleh Dindar and the Swarm-NG Development Team  *
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

/*! 
 *  \file hermite_bpt.hpp Body-per-thread implementation of Hermite PEC2 integrator
 */

#include "swarm/common.hpp"
#include "swarm/gpu/bppt.hpp"
#include "vector_types.h"
#include "vector_functions.h"

namespace swarm { namespace gpu {


GENERIC double3 operator+(const double3 &a,const double3 &b) {
  return make_double3(a.x+b.x, a.y+b.y, a.z+b.z);
}
GENERIC double3 operator-(const double3 &a,const double3 &b) {
  return make_double3(a.x-b.x, a.y-b.y, a.z-b.z);
}
GENERIC double3 operator*(const double3 &a, const double& b) {
  return make_double3(a.x*b, a.y*b, a.z*b);
}
GENERIC double3 operator*(const double& b, const double3& a) {
  return a * b;
}
GENERIC double3 operator/(const double3 &a, const double& b){
  return a *(1.0/b);
}

GENERIC double3& operator+=(double3& a, const double3& b){
  a = a + b;
}
GENERIC double3& operator-=(double3& a, const double3& b){
  a = a - b;
}

GENERIC double sqnorm(double3 a) { return a.x*a.x+a.y*a.y+a.z*a.z; }
GENERIC double inner_product(double3 a,double3 b){ return a.x*b.x+a.y*b.y+a.z*b.z; }

template<class T>
class BodyGravitation {
  static const int nbod = T::n;
  ensemble::SystemRef& sys;
  const int& b;
public:
  GENERIC BodyGravitation(ensemble::SystemRef& sys,const int& b):sys(sys),b(b){}
  
  GENERIC void operator()(const double3& pos, const double3& vel, double3& acc,double3& jerk) {
    // Look over all bodies saving sun for last
    acc = make_double3(0.0,0.0,0.0);
    jerk = make_double3(0.0,0.0,0.0);
    
#pragma unroll
    for(int d = nbod-1; d >= 0; d--) 
      if( d != b ) {
        double M  = sys[d].mass();
        double3 dx = sys[d].pos() - pos;
        double3 dv = sys[d].vel() - vel;
        double r2 = sqnorm(dx), rinv = rsqrt(r2) / r2, rv = inner_product(dx,dv) * 3.0 / r2;
        acc +=  M * rinv * dx;
        jerk += M * rinv * ( dv - rv * dx );
      }
  }
  
};

  
/*! A self-contained GPU integrator that is independent of
 * most of the framework. Just to test how much advantage do 
 * we get by using body-per-thread distribution of the workload
 * \ingroup integrators
 */
template<class Monitor, template<class T> class Gravitation>
class hermite_bpt : public integrator {
  typename Monitor::params _monitor_params;
  double _time_step;
  int _spb;
public:
  explicit hermite_bpt(const config& cfg):integrator(cfg),_monitor_params(cfg) {
    _time_step = cfg.require("time_step", 0.001);
    _spb = cfg.optional("system_per_block", 0);
  }
  
  virtual void launch_integrator() {
    launch_templatized_integrator(this);
  }
  
  template<class T>
  static GENERIC int thread_per_system(T ctp) {
    return T::n;
  }
  
  template<class T>
  static GENERIC int shmem_per_system(T ctp) {
    return 0;
  }
  
  int override_system_per_block(){
    return _spb;
  }
  
  
  template<class T>
  GPUAPI void kernel(T ctp){
    if(sysid() >= _dens.nsys()) return;
    
    ensemble::SystemRef sys = _dens[sysid()];
    const int nbod = T::n;
    const int b = thread_in_system();
    Gravitation<T> calcForces(sys, b);
    
    //Monitor monitor_test(_monitor_params,sys,*_log);
    
    double3 pos = sys[b].pos(), vel = sys[b].vel();
    
    
    double3 acc0, jerk0;
    calcForces(pos, vel, acc0, jerk0);
    
    for(int iter = 0; (iter < _max_iterations) && sys.is_active(); iter++ ) {

      double h;
      if( sys.time() + _time_step < _destination_time )
        h = _time_step;
      else
        h = _destination_time - sys.time();
      
      // Predict
      pos = pos + h * (vel + (h * 0.5) * (acc0 + (h / 3.0) * jerk0));
      vel = vel + h * (acc0 + (h * 0.5) * jerk0);
      
      
      double3 acc1, jerk1;
      {
        // Evaluate
        sys[b].set_pos(pos), sys[b].set_vel(vel);__syncthreads();
        calcForces(pos, vel, acc1, jerk1);
        
        // Correct
        pos += ( (0.1-0.25) * (acc0 - acc1) - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h) * h * h;
        vel += (( -0.5 ) * (acc0 - acc1 ) -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h )* h ;
      }
      {
        // Evaluate
        sys[b].set_pos(pos), sys[b].set_vel(vel);__syncthreads();
        calcForces(pos, vel, acc1, jerk1);
        
        // Correct
        pos += ( (0.1-0.25) * (acc0 - acc1) - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h) * h * h;
        vel += (( -0.5 ) * (acc0 - acc1 ) -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h )* h ;
      }
      acc0 = acc1, jerk0 = jerk1;
      
      // Store the position and velocity back to global memory
      sys[b].set_pos(pos), sys[b].set_vel(vel);
      
      if( b == 0 ) {
        sys.time() += h;
        if( sys.time() >= _destination_time ) sys.set_inactive();
      }
      __syncthreads();
      
      //monitor_test();
      
    }
  }
};
  
} }