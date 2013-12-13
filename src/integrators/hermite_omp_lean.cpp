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
            double r2 = sqnorm(dx), rinv = 1/sqrt(r2) / r2, rv = inner_product(dx,dv) * 3.0 / r2;
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

} // namespace swarm


