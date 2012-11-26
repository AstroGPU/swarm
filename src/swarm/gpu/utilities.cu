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

/*! \file utilities.cu
 *  \brief provides several utility  functions for public interface for swarm libaray
 *
 */

#include "../common.hpp"
#include "utilities.hpp"


namespace swarm {

struct count_systems_t {
	deviceEnsemble ens;
	int count_running;

};

__global__ 
void number_of_active_systems_kernel( count_systems_t* csys){
	int sysid = ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	if(sysid>=csys->ens.nsys()) return;
	if(csys->ens[sysid].is_active() )
		atomicAdd(&csys->count_running,1);
};

int number_of_active_systems(deviceEnsemble ens) {
	const int system_per_block = 16*16;  // need not be the same as used by the integrator
	const int nblocks = ( ens.nsys() + system_per_block - 1 ) / system_per_block;
	dim3 gD; gD.z = 1;
	find_best_factorization(gD.x,gD.y,nblocks);
	dim3 tD; tD.x = system_per_block; tD.y = 1;

	count_systems_t count_systems, *pcount_systems ;

	count_systems.ens = ens;
	count_systems.count_running = 0;


	cudaErrCheck ( cudaMalloc(&pcount_systems,sizeof(count_systems_t)) );
	cudaErrCheck ( cudaMemcpy(pcount_systems,&count_systems,sizeof(count_systems_t),cudaMemcpyHostToDevice) );

	number_of_active_systems_kernel<<< gD, tD >>>( pcount_systems );

	cudaErrCheck ( cudaMemcpy(&count_systems,pcount_systems,sizeof(count_systems_t),cudaMemcpyDeviceToHost) );
	cudaErrCheck ( cudaFree(pcount_systems) );

	return count_systems.count_running;

}

bool configure_grid(dim3 &gridDim, int threadsPerBlock, int nthreads, int dynShmemPerThread, int staticShmemPerBlock)
{
        const int shmemPerMP =  16384;

        int dyn_shared_mem_required = dynShmemPerThread*threadsPerBlock;
        int shared_mem_required = staticShmemPerBlock + dyn_shared_mem_required;
        if(shared_mem_required > shmemPerMP) { return false; }

        // calculate the total number of threads
        int nthreadsEx = nthreads;
        int over = nthreads % threadsPerBlock;
        if(over) { nthreadsEx += threadsPerBlock - over; } // round up to multiple of threadsPerBlock

        // calculate the number of blocks
        int nblocks = nthreadsEx / threadsPerBlock;
        if(nthreadsEx % threadsPerBlock) { nblocks++; }

        // calculate block dimensions so that there are as close to nblocks blocks as possible
        find_best_factorization(gridDim.x, gridDim.y, nblocks);
        gridDim.z = 1;

        return true;
}



void find_best_factorization(unsigned int &bx, unsigned int &by, int nblocks)
{
        bx = -1;
        int best_r = 100000;
        for(int bytmp = 1; bytmp != 65536; bytmp++)
        {
                int r  = nblocks % bytmp;
                if(r < best_r && nblocks / bytmp < 65535)
                {
                        by = bytmp;
                        bx = nblocks / bytmp;
                        best_r = r;

                        if(r == 0) { break; }
                        bx++;
                }
        }
        if(bx == -1) { std::cerr << "Unfactorizable?!\n"; exit(-1); }
}


__global__ 
void reactivate_systems_kernel( deviceEnsemble ens ){
	int sysid = ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	if(sysid>=ens.nsys()) return;
	if(ens[sysid].is_inactive() )
	   ens[sysid].set_active();
};



void reactivate_systems(deviceEnsemble ens) {
	const int system_per_block = 16*16;  // need not be the same as used by the integrator
	const int nblocks = ( ens.nsys() + system_per_block - 1 ) / system_per_block;
	dim3 gD; gD.z = 1;
	find_best_factorization(gD.x,gD.y,nblocks);
	dim3 tD; tD.x = system_per_block; tD.y = 1;

	reactivate_systems_kernel<<< gD, tD >>>( ens );
}


}
