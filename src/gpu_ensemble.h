/*************************************************************************
 * Copyright (C) 2010 by Mario Juric  and the Swarm-NG Development Team  *
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

/*! \file gpu_ensemble.h
 *   \brief Declares gpu_ensemble class
 *
*/

/// The main namespace for the Swarm-NG library
namespace swarm {
/**
        \brief Specialization of ensemble, with data stored on the device (GPU)

        Use this class to create/read/write ensembles in GPU memory. Typically, you
        will want to create a cpu_ensemble, populate it on the host, use gpu_ensemble's
        copy constructor to upload the data to GPU, and pass the gpu_ensemble structure
        to an integrator that will advance the ensemble.
*/
class gpu_ensemble : public ensemble
{
        protected:
                int *nactive_gpu;       //!< temp variable for get_nactive()

        public:
		//constructors
                gpu_ensemble();                                 //!< instantiate an empty ensemble
                gpu_ensemble(int nsys, int nbod);               //!< instantiate an ensemble with room for nsys systems with nbod each
                gpu_ensemble(const cpu_ensemble &source);       //!< instantiate a copy of the source ensemble

                ~gpu_ensemble();

		//Methods defined in gpu_ensemble.cpp
                void copy_from(const cpu_ensemble &source);     //!< Copy the data from the CPU
                void free();                                                    //!< Deallocate GPU memory
                void reset(int nsys, int nbod, bool reinitIndices = true);      //!< Allocate GPU memory for nsys systems of nbod planets each

		//Methods defined in gpu_ensemble.cu 
		/* TODO: All of these methods are currently commented out of gpu_ensemble.cu */
                int get_nactive() const;                        //!< Download and return the number of active systems
                void set_time_end_all(const real_time tend);
                __device__ void set_time_end_all_kernel(const real_time tend);
};
//end class gpu_ensemble
}
