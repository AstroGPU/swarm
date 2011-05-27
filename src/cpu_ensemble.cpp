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

/*! \file cpu_ensemble.cpp
 *  \brief class cpu_ensemble 
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
   \brief CPU Ensemble class plumbing (mostly memory management)
*/
cpu_ensemble::cpu_ensemble()
{
        construct_base();
}

/*!
   \brief CPU Ensemble class plumbing (mostly memory management)

  @param[in] sys number of systems
  @param[in] nbod number of bodies
*/
cpu_ensemble::cpu_ensemble(int nsys, int nbod)
{
        construct_base();
        reset(nsys, nbod);
}

/*!
   \brief CPU Ensemble class plumbing (mostly memory management)

  @param[in] source gpu_ensemble
*/
cpu_ensemble::cpu_ensemble(const gpu_ensemble &source)
{
        construct_base();
        copy_from(source);
}

/*!
   \brief CPU Ensemble class plumbing (mostly memory management)

  @param[in] source cpu_ensemble
*/
cpu_ensemble::cpu_ensemble(const cpu_ensemble &source)
{
        construct_base();
        copy_from(source);
}


/********** METHODS DEFINED IN cpu_ensemble.h **********/

/*!
   \brief Copy the data from the CPU

  @param[in] src cpu_ensemble
*/
void cpu_ensemble::copy_from(const cpu_ensemble &src)
{
        reset(src.nsys(), src.nbod(), false);

        // low-level copy from host to separate host memory
        memcpy(m_T, src.m_T, m_nsys*sizeof(*m_T));
        memcpy(m_Tend, src.m_Tend, m_nsys*sizeof(*m_Tend));
        memcpy(m_nstep, src.m_nstep, m_nsys*sizeof(*m_nstep));
        memcpy(m_Toutput, src.m_Toutput, 2*m_nsys*sizeof(*m_Toutput));
        memcpy(m_xyz, src.m_xyz, 3*m_nbod*m_nsys*sizeof(*m_xyz));
        memcpy(m_vxyz, src.m_vxyz, 3*m_nbod*m_nsys*sizeof(*m_vxyz));
        memcpy(m_m, src.m_m, m_nbod*m_nsys*sizeof(*m_m));
        memcpy(m_flags, src.m_flags, m_nsys*sizeof(*m_flags));
        memcpy(m_systemIndices, src.m_systemIndices, m_nsys*sizeof(*m_systemIndices));
//      memcpy(m_nactive, src.m_nactive, 1*sizeof(*m_nactive));

        m_last_integrator = src.m_last_integrator;
}


/*!
   \brief Copy the data from the GPU

  @param[in] src gpu_ensemble
*/
void cpu_ensemble::copy_from(const gpu_ensemble &src)
{
        reset(src.nsys(), src.nbod(), false);

        // low-level copy from host to device memory
        memcpyToHost(m_T, src.m_T, m_nsys);
        memcpyToHost(m_Tend, src.m_Tend, m_nsys);
        memcpyToHost(m_nstep, src.m_nstep, m_nsys);
        memcpyToHost(m_Toutput, src.m_Toutput, 2*m_nsys);
        memcpyToHost(m_xyz, src.m_xyz, 3*m_nbod*m_nsys);
        memcpyToHost(m_vxyz, src.m_vxyz, 3*m_nbod*m_nsys);
        memcpyToHost(m_m, src.m_m, m_nbod*m_nsys);
        memcpyToHost(m_flags, src.m_flags, m_nsys);
        memcpyToHost(m_systemIndices, src.m_systemIndices, m_nsys);
//      memcpyToHost(m_nactive, src.m_nactive, 1);

        m_last_integrator = src.m_last_integrator;
}


/*!
   \brief  Deallocate CPU memory
*/
void cpu_ensemble::free()
{
//      hostFree(m_nactive); m_nactive = NULL;
        hostFree(m_T); m_T = NULL;
        hostFree(m_Tend); m_Tend = NULL;
        hostFree(m_nstep); m_nstep = NULL;
        hostFree(m_Toutput); m_Toutput = NULL;
        hostFree(m_xyz); m_xyz = NULL;
        hostFree(m_vxyz); m_vxyz = NULL;
        hostFree(m_m); m_m = NULL;
        hostFree(m_flags); m_flags = NULL;
        hostFree(m_systemIndices); m_systemIndices = NULL;
}

/*!
   \brief data packing for CPU ensemble

   ...
  @param[in] src cpu_ensemble
  @return ...
*/
int cpu_ensemble::pack()
{
    int openid=0;
    for(int sysid=0;sysid<nsys();++sysid)
      {
        if(is_active(sysid))
          {
            if(sysid!=openid)
              {
                std::swap(m_T[sysid],m_T[openid]);
                std::swap(m_Tend[sysid],m_Tend[openid]);
                std::swap(m_nstep[sysid],m_nstep[openid]);
                std::swap(m_Toutput[sysid],m_Toutput[openid]);
                std::swap(m_Toutput[sysid+nsys()],m_Toutput[openid+nsys()]);
                std::swap(m_flags[sysid],m_flags[openid]);
                std::swap(m_systemIndices[sysid],m_systemIndices[openid]);
                for(int bod=0;bod<nbod();++bod)
                  {
                    double tmp;
                    tmp = mass(openid,bod); mass(openid,bod) = mass(sysid,bod); mass(sysid,bod) = tmp;
                    tmp = x(openid,bod);   x(openid,bod) =  x(sysid,bod);  x(sysid,bod) = tmp;
                    tmp = y(openid,bod);   y(openid,bod) =  y(sysid,bod);  y(sysid,bod) = tmp;
                    tmp = z(openid,bod);   z(openid,bod) =  z(sysid,bod);  z(sysid,bod) = tmp;
                    tmp = vx(openid,bod); vx(openid,bod) = vx(sysid,bod); vx(sysid,bod) = tmp;
                    tmp = vy(openid,bod); vy(openid,bod) = vy(sysid,bod); vy(sysid,bod) = tmp;
                    tmp = vz(openid,bod); vz(openid,bod) = vz(sysid,bod); vz(sysid,bod) = tmp;
                  }
              } // end if sysid!=openid
                          ++openid;
          } // end if is_active
      }  // end for sysid
    return openid;
}

/*!
   \brief Merge in systems from another cpu_ensemble

   ...
  @param[in,out] src cpu_ensemble  ...
  @param[in] offset ...
*/
void cpu_ensemble::replace_inactive_from(cpu_ensemble &src, const int offset)   // Merge in systems from another cpu_ensemble
{
        int src_sys = 0, dest_sys = 0;
        //      for(int dest_sys=start_sys;dest_sys<nsys();++dest_sys)
        do
          {
            //      do while(is_inactive(dest_sys)&&(dest_sys<nsys())) { ++dest_sys; };
            for(;is_inactive(dest_sys)&&(dest_sys<nsys());++dest_sys);
            if(dest_sys>=nsys()) break;
            do { ++src_sys; } while (src.is_inactive(src_sys)&&(src_sys<src.nsys()));
            if(src_sys>=src.nsys()) break;
            assert(src.is_active(src_sys));
            time(dest_sys) = src.time(src_sys);
            time_end(dest_sys) = src.time_end(src_sys);
            nstep(dest_sys) = src.nstep(src_sys);
            time_output(dest_sys,0) = src.time_output(src_sys,0);
            time_output(dest_sys,1) = src.time_output(src_sys,1);
            flags(dest_sys) = src.flags(src_sys);
            index_of_system(dest_sys) = src.index_of_system(src_sys) + offset;
            for(int bod=0;bod<nbod();++bod)
              {
                mass(dest_sys,bod) = mass(src_sys,bod);
                x (dest_sys,bod) = x (src_sys,bod);
                y (dest_sys,bod) = y (src_sys,bod);
                z (dest_sys,bod) = z (src_sys,bod);
                vx(dest_sys,bod) = vx(src_sys,bod);
                vy(dest_sys,bod) = vy(src_sys,bod);
                vz(dest_sys,bod) = vz(src_sys,bod);
              }
            src.set_inactive(src_sys);
          } while(true);  // end loop over dest_sys

        m_last_integrator = NULL;  // I beleive this is necessary since we don't know that other systems have been integrated the same way.  Any problem?
}


/*!
   \brief  Allocate CPU memory for nsys systems of nbod planets each

  @param[in] sys number of systems
  @param[in] nbod number of bodies
  @param[in] reinitIndices flag for reinitialize indices
*/
void cpu_ensemble::reset(int nsys, int nbod, bool reinitIndices)        //
{
        // do we need to realloc?
        if(m_nsys != nsys || m_nbod != nbod)
        {
                m_T = hostAlloc(m_T, nsys);
                m_Tend = hostAlloc(m_Tend, nsys);
                m_nstep = hostAlloc(m_nstep, nsys);
                m_Toutput = hostAlloc(m_Toutput, 2*nsys);
                m_xyz = hostAlloc(m_xyz, 3*nsys*nbod);
                m_vxyz = hostAlloc(m_vxyz, 3*nsys*nbod);
                m_m = hostAlloc(m_m, nsys*nbod);
                m_flags = hostAlloc(m_flags, nsys);
                m_systemIndices = hostAlloc(m_systemIndices, nsys);

//              m_nactive = hostAlloc(m_nactive, 1);
                m_nsys = nsys;
                m_nbod = nbod;
        }

        if(reinitIndices)
        {
                // reinitialize system indices
                for(int i = 0; i != nsys; i++) { m_systemIndices[i] = i; }

                // Set all systems to active
                for(int sys = 0; sys != nsys; sys++) { flags(sys) = 0; }

                // Set nstep=0
                for(int sys = 0; sys != nsys; sys++) { nstep(sys) = 0; }
        }

        m_last_integrator = NULL;
}

}
