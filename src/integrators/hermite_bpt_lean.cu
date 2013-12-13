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
#include "double3.hpp"



namespace swarm { 



namespace gpu {

template<class T>
__global__ void kernel_launch(T self){
    return self.kernel();
}

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

//        std::cout << gridDim.x << " " << gridDim.y << " " << gridDim.z << " "
//            << " " << threadDim.x << " " << threadDim.y << " " << threadDim.z << std::endl;

        kernel_launch<<<gridDim, threadDim, 0>>>(*this);
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
