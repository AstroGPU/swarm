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

/*! \file gpu_ensemble.cu 
 *  \brief cuda code for gpu_kernel. 
 *  TODO: This ic currently all commented out by a preprocessor directive!
*/

#include "swarm.h"
#include <cux/cux.h>
#include <astro/util.h>
#include <cassert>

namespace swarm {

#if 0

TODO: Need to dynamically compute grid size for this to work properly.
      Do not enable until this is done.
int gpu_ensemble::get_nactive() const
{
        if(nactive_gpu)
        {
                cudaMalloc((void**)&nactive_gpu, sizeof(*nactive_gpu));
        }

        cudaMemset(nactive_gpu, 0, sizeof(*nactive_gpu));
        get_nactive_kernel<<<60, 64>>>(nactive_gpu, *this);

        // fetch the result
        int nrunning;
        memcpyToHost(&nrunning, nactive_gpu, 1);
        return nrunning;
}
#endif

#if 0  // compiler says not allowed, need to figure out
__device__ void gpu_ensemble::set_time_end_all_kernel(const real_time tend)
{
int sys = threadId();
if(sys>= nsys()) { return; }
time_end(sys) = tend;
}

void gpu_ensemble::set_time_end_all(const real_time tend)
{
  set_time_end_all_kernel<<<60, 64>>>(tend);
}
#endif
}
