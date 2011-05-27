/*************************************************************************
 * Copyright (C) 2008-2010 by Mario Juric & Swarm-NG Development Team    *
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

/*! \file gpu_ensemble.cpp
 *  \brief class gpu_ensemble
 *
*/

#include <cuda_runtime_api.h>
#include <algorithm> // for swap
#include <memory>
#include <iostream>
#include <dlfcn.h>
#include <fstream>
#include "swarm.h"

namespace swarm {

/*************** CONSTRUCTORS ****************/
/*!
   \brief GPU Ensemble class constructor
*/
gpu_ensemble::gpu_ensemble()
{
        construct_base();
        nactive_gpu = NULL;
}

/*!
   \brief GPU Ensemble class constructor

  @param[in] sys number of systems
  @param[in] nbod number of bodies
*/
gpu_ensemble::gpu_ensemble(int nsys, int nbod)
{
        construct_base();
        reset(nsys, nbod);
        nactive_gpu = NULL;
}

/*!
   \brief GPU Ensemble class constructor

  @param[in] source cpu_ensemble
*/
gpu_ensemble::gpu_ensemble(const cpu_ensemble &source)
{
        construct_base();
        copy_from(source);
        nactive_gpu = NULL;
}

/*!
   \brief GPU Ensemble class destructor

   Deallocate GPU memory.
*/
gpu_ensemble::~gpu_ensemble()
{
        free();
        cudaFree(nactive_gpu);
}


/********** METHODS DEFINED IN gpu_ensemble.h **********/

/*!
   \brief Copy the data from the CPU

  @param[in] src cpu_ensemble
*/
void gpu_ensemble::copy_from(const cpu_ensemble &src)
{
        reset(src.nsys(), src.nbod(), false);

        // low-level copy from host to device memory
        memcpyToGPU(m_T, src.m_T, m_nsys);
        memcpyToGPU(m_Tend, src.m_Tend, m_nsys);
        memcpyToGPU(m_nstep, src.m_nstep, m_nsys);
        memcpyToGPU(m_Toutput, src.m_Toutput, 2*m_nsys);
        memcpyToGPU(m_xyz, src.m_xyz, 3*m_nbod*m_nsys);
        memcpyToGPU(m_vxyz, src.m_vxyz, 3*m_nbod*m_nsys);
        memcpyToGPU(m_m, src.m_m, m_nbod*m_nsys);
        memcpyToGPU(m_flags, src.m_flags, m_nsys);
        memcpyToGPU(m_systemIndices, src.m_systemIndices, m_nsys);
//      memcpyToGPU(m_nactive, src.m_nactive, 1);

        m_last_integrator = src.m_last_integrator;
}

/*!
   \brief  Deallocate GPU memory
*/
void gpu_ensemble::free()
{
//      cudaFree(m_nactive); m_nactive = NULL;
        cudaFree(m_T); m_T = NULL;
        cudaFree(m_Tend); m_Tend = NULL;
        cudaFree(m_nstep); m_nstep = NULL;
        cudaFree(m_Toutput); m_Toutput = NULL;
        cudaFree(m_xyz); m_xyz = NULL;
        cudaFree(m_vxyz); m_vxyz = NULL;
        cudaFree(m_m); m_m = NULL;
        cudaFree(m_flags); m_flags = NULL;
        cudaFree(m_systemIndices); m_systemIndices = NULL;
}

/*!
   \brief  Allocate GPU memory for nsys systems of nbod planets each

  @param[in] sys number of systems
  @param[in] nbod number of bodies
  @param[in] reinitIndices flag for reinitialize indices
*/
void gpu_ensemble::reset(int nsys, int nbod, bool reinitIndices)
{
        // do we need to realloc?
        if(m_nsys != nsys || m_nbod != nbod)
        {
                free();

//              cudaMalloc((void**)&m_nactive, sizeof(*m_nactive));
                cudaMalloc((void**)&m_T, nsys*sizeof(*m_T));
                cudaMalloc((void**)&m_Tend, nsys*sizeof(*m_Tend));
                cudaMalloc((void**)&m_nstep, nsys*sizeof(*m_nstep));
                cudaMalloc((void**)&m_Toutput, 2*nsys*sizeof(*m_Toutput));
                cudaMalloc((void**)&m_xyz, 3*nsys*nbod*sizeof(*m_xyz));
                cudaMalloc((void**)&m_vxyz, 3*nsys*nbod*sizeof(*m_vxyz));
                cudaMalloc((void**)&m_m, nsys*nbod*sizeof(*m_m));
                cudaMalloc((void**)&m_flags, nsys*sizeof(*m_flags));
                cudaMalloc((void**)&m_systemIndices, nsys*sizeof(*m_systemIndices));

                m_nsys = nsys;
                m_nbod = nbod;
        }

        if(reinitIndices)
        {
                // reinitialize system indices
                std::vector<int> tmp(nsys);
                for(int i = 0; i != nsys; i++) { tmp[i] = i; }
                cudaMemcpy(m_systemIndices, &tmp[0], tmp.size()*sizeof(tmp[0]), cudaMemcpyHostToDevice);

                // Set all systems active
                std::valarray<int> tmp2(nsys);
                tmp2 = 0;
                cudaMemcpy(m_flags, &tmp2[0], tmp2.size()*sizeof(tmp2[0]), cudaMemcpyHostToDevice);

                // Set nsteps = 0
                std::valarray<uint> tmp3(nsys);
                tmp3 = 0;
                cudaMemcpy(m_nstep, &tmp3[0], tmp3.size()*sizeof(tmp3[0]), cudaMemcpyHostToDevice);

                // nactive
//              int tmp4 = nsys;
//              cudaMemcpy(m_nactive, &tmp4, sizeof(tmp4), cudaMemcpyHostToDevice);
        }

        // clear the m_last_integrator field
        m_last_integrator = NULL;
}

}
