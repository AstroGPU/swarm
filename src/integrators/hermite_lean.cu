/*************************************************************************
 * Copyright (C) 2013 by Saleh Dindar and the Swarm-NG Development Team  *
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
/*! \file hermite_lean.cpp
 *  \brief A lean and simple multi-threaded Hermite PEC2 integrator on CPU and GPU
 *
 *  This file is a simplified version of \ref hermite_cpu.hpp and 
 *  \ref hermite_bpt.hpp just for 
 *  updated benchmark. This is also a basis for furthur refactoring of
 *  the rest of the code
 */

#include "swarm/common.hpp"
#include "swarm/integrator.hpp"
#include "swarm/plugin.hpp"
#include "swarm/gpu/device_settings.hpp"
#include "vector_types.h"
#include "vector_functions.h"

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


namespace swarm { 



namespace cpu {

/*! Multi-threaded implementation of Hermite PEC2
 *
 */
class hermite_omp_lean : public integrator {
    double _time_step;

public:
    hermite_omp_lean(const config& cfg):integrator(cfg) {
        _time_step = cfg.require("time_step", 0.0);
    }

    virtual void launch_integrator(){
        #pragma omp parallel for
        for(int i = 0; i < _ens.nsys(); i++)
            integrate_system(_ens[i]);
    }


    void calcForces(ensemble::SystemRef& sys, double3 acc[], double3 jerk[]){
        const int nbod = sys.nbod();

        for(int i = 0; i < nbod; i++)
            acc[i] = make_double3(0.0,0.0,0.0), jerk[i] = make_double3(0.0,0.0,0.0);

        for(int i = 0; i < nbod-1; i++) for(int j = i+1; j < nbod; j++) {
            double3 dx = sys[j].pos() - sys[i].pos(),
                    dv = sys[j].vel() - sys[i].vel();
            double r2 = sqnorm(dx), rinv = rsqrt(r2) / r2, rv = inner_product(dx,dv) * 3.0 / r2;
            acc [i] += sys[j].mass() * rinv * dx;
            jerk[i] += sys[j].mass() * rinv * ( dv - rv * dx );

            acc [j] -= sys[i].mass() * rinv * dx;
            jerk[j] -= sys[i].mass() * rinv * ( dv - rv * dx );
        }
    }

    void integrate_system(ensemble::SystemRef sys) {
        const int nbod = sys.nbod();

        double3 acc0[nbod], acc1[nbod], jerk0[nbod], jerk1[nbod];
        calcForces(sys, acc0, jerk0);


		for(int iter = 0 ; (iter < _max_iterations) && sys.is_active() ; iter ++ ) {
			double h = _time_step;

			if( sys.time() + h > _destination_time ) {
				h = _destination_time - sys.time();
			}

            // Predict
            for(int b = 0; b < nbod; b++) {
                sys[b].set_pos(sys[b].pos() + h * (sys[b].vel() + (h * 0.5) * (acc0[b] + (h/3.0)* jerk0[b])));
                sys[b].set_vel(sys[b].vel() + h * (acc0[b] + h * 0.5 * jerk0[b]));
            }


            // Evaluate and Correct round 1
            {
                calcForces(sys,acc1,jerk1);
                for(int b = 0; b < nbod; b++) {
                    sys[b].set_pos(
                    sys[b].pos() + (.1-.25)* (acc0[b] - acc1[b]) * h * h
                        - 1/60.0 * ( 7 * jerk0[b] + 2 * jerk1[b]) * h * h * h
                        );
                    sys[b].set_vel(
                    sys[b].vel() - 0.5 * (acc0[b] - acc1[b]) * h
                        - 1/12.0 * ( 5 * jerk0[b] + jerk1[b] ) * h * h
                        );
                }
            }
            // Evaluate and Correct round 2
            {
                calcForces(sys,acc1,jerk1);
                for(int b = 0; b < nbod; b++) {
                    sys[b].set_pos(
                    sys[b].pos() + (.1-.25)* (acc0[b] - acc1[b]) * h * h
                        - 1/60.0 * ( 7 * jerk0[b] + 2 * jerk1[b]) * h * h * h
                        );
                    sys[b].set_vel(
                    sys[b].vel() - 0.5 * (acc0[b] - acc1[b]) * h
                        - 1/12.0 * ( 5 * jerk0[b] + jerk1[b] ) * h * h
                        );
                }
            }

            for(int b = 0; b < nbod; b++)
                acc0[b] = acc1[b], jerk0[b] = jerk1[b];

            sys.time() += h;

            if( sys.time() > _destination_time - 1e-12 )
                sys.set_inactive();

        }
    }
    
};

integrator_plugin_initializer<hermite_omp_lean> hermite_omp_lean_plugin("hermite_omp_lean");

}  // namespace cpu


template<class T>
__global__ void kernel_launch(T* self){
    return self->kernel();
}


namespace gpu {

class hermite_bpt_lean : public integrator {
    double _time_step;
    //! Number of systems per block
    int _spb;

public:
    explicit hermite_bpt_lean(const config& cfg):integrator(cfg) {
        _time_step = cfg.require("time_step", 0.001);
        _spb = cfg.optional("system_per_block", 0);
    }

    virtual void launch_integrator(){
        int bb = 1;
        const int spb = _spb > 0 ? _spb
            : optimized_system_per_block(SHMEM_CHUNK_SIZE, _ens.nbod(), 0, bb);
        
        dim3 threadDim;
        threadDim.x = spb;
        threadDim.y = _ens.nbod();
        threadDim.z = 1;

		if(!check_cuda_limits(threadDim.x * threadDim.y * threadDim.z, 0))
			throw runtime_error("The block size settings exceed CUDA requirements");

        const int nblocks = (_ens.nsys() + spb - 1 ) / spb;
        dim3 gridDim;
        gridDim.x = nblocks; gridDim.y = 1; gridDim.z = 1;

        std::cout << gridDim.x << " " << gridDim.y << " " << gridDim.z << " "
            << " " << threadDim.x << " " << threadDim.y << " " << threadDim.z << std::endl;

        typedef hermite_bpt_lean implementation;
        implementation *integ = this, *gpu_integ;
		cudaErrCheck ( cudaMalloc(&gpu_integ,sizeof(implementation)) );
		cudaErrCheck ( cudaMemcpy(gpu_integ,integ,sizeof(implementation),cudaMemcpyHostToDevice) );

        kernel_launch<<<gridDim, threadDim, 0>>>(gpu_integ);

        cudaErrCheck ( cudaFree(gpu_integ) );
    }


    static GPUAPI void calcForces(ensemble::SystemRef& sys, const int b, double3& acc, double3& jerk) {

        const int nbod = sys.nbod();
        acc = make_double3(0.0,0.0,0.0), jerk = make_double3(0.0,0.0,0.0);
        for(int d = nbod-1; d >= 0; d--) if(d != b) {
            double M  = sys[d].mass();
            double3 dx = sys[d].pos() - sys[b].pos();
            double3 dv = sys[d].vel() - sys[b].vel();
            double r2 = sqnorm(dx), rinv = rsqrt(r2) / r2, rv = inner_product(dx,dv) * 3.0 / r2;
            acc +=  M * rinv * dx;
            jerk += M * rinv * ( dv - rv * dx );
        }
    }


    GPUAPI void kernel() {
      const int blockNumber = (blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x;
      const int sysid = (blockNumber * blockDim.z + threadIdx.z) * blockDim.x + threadIdx.x;
      const int b = threadIdx.y;

      if(sysid >= _dens.nsys()) return;
      ensemble::SystemRef sys = _dens[sysid];
      const int nbod = sys.nbod();

      double3 pos = sys[b].pos(), vel = sys[b].vel();

      double3 acc0, jerk0, acc1, jerk1;
      calcForces(sys, b, acc0, jerk0);

        for(int iter = 0; (iter < _max_iterations) && sys.is_active(); iter++ ) {

          double h;
          if( sys.time() + _time_step < _destination_time )
            h = _time_step;
          else
            h = _destination_time - sys.time();

          // Predict
          pos = pos + h * (vel + (h * 0.5) * (acc0 + (h / 3.0) * jerk0));
          vel = vel + h * (acc0 + (h * 0.5) * jerk0);

          //// Round 1
          {
            // Evaluate
            sys[b].set_pos(pos), sys[b].set_vel(vel);__syncthreads();
            calcForces(sys, b,acc1, jerk1);
            
            // Correct
            pos += ( (0.1-0.25) * (acc0 - acc1) - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h) * h * h;
            vel += (( -0.5 ) * (acc0 - acc1 ) -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h )* h ;
          }
          //// Round 2
          {
            // Evaluate
            sys[b].set_pos(pos), sys[b].set_vel(vel);__syncthreads();
            calcForces(sys, b, acc1, jerk1);
            
            // Correct
            pos += ( (0.1-0.25) * (acc0 - acc1) - 1.0/60.0 * ( 7.0 * jerk0 + 2.0 * jerk1 ) * h) * h * h;
            vel += (( -0.5 ) * (acc0 - acc1 ) -  1.0/12.0 * ( 5.0 * jerk0 + jerk1 ) * h )* h ;
          }

          acc0 = acc1, jerk0 = jerk1;

          // Store position and velocity to memory
          sys[b].set_pos(pos), sys[b].set_vel(vel);

          // The zeroth thread should update the time
          if( b == 0 ) {
            sys.time() += h;
            if( sys.time() >= _destination_time ) sys.set_inactive();
          }
          __syncthreads();


        }

    }
};

integrator_plugin_initializer<hermite_bpt_lean> hermite_bpt_lean_plugin("hermite_bpt_lean");


} // namespace gpu



} // namespace swarm
