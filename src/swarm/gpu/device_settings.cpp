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

/*! \file device_settings.cpp
 *   \brief Implements the function to set up the GPU related
 *          parameters. 
 *
 */

#include "swarm/common.hpp"
#include "device_settings.hpp"

#if 0
// Setting prefered cache size the following function
// can do this but the problem is that since our kernel is really
// a template functin we don't know the exact name for it
void set_prefered_shared_memory(const char* function_name){
	cudaFuncSetCacheConfig( name, cudaFuncCachePreferShared);
}
#endif

const int registers_per_thread = 64;  
cudaDeviceProp deviceInfo;

/*
 * \todo Is it intentional that shmem_per_system isn't multiplied by chunk_size?
 */
int optimized_system_per_block(int chunk_size, int thread_per_system
		, int shmem_per_system){
	return blocks_per_mp( chunk_size * thread_per_system, shmem_per_system) 
		* chunk_size ;
}

void select_cuda_device(int dev) {
	int devcnt; cudaErrCheck( cudaGetDeviceCount(&devcnt) );
	if( dev >= 0 && dev < devcnt )
		cudaErrCheck( cudaSetDevice(dev) );
	else
		std::cerr << "Cannot select the CUDA device. GPU integrators are disabled" << std::endl;

	  cudaErrCheck( cudaGetDeviceProperties(&deviceInfo, dev) );


}


void set_more_cache(){
	$$$;
	cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
}

void print_device_information(){
	  std::cerr << "Device:\t"  << deviceInfo.name   << "\n"
		  //<< " Compute Capabality: " << deviceInfo.computeMode <<   "\n"
		  << "Global Memory:\t" << deviceInfo.totalGlobalMem/double(1<<30)     << "GB\n"
		  << "Shared Memory\t"  << deviceInfo.sharedMemPerBlock/1024  << "KB\n"
		  << "Max Blocksize\t"  << deviceInfo.maxThreadsPerBlock << "\n"
		  << "Warp Size    \t"  << deviceInfo.warpSize   << "\n"
		  << "Registers/MP \t"  << deviceInfo.regsPerBlock       << "\n"
		  ;
	  
}

/*
 * \todo is block_warps computed correctly when blocksize is a multiple of warpSize?
 */
int blocks_per_mp( int blocksize, int shmem_per_block ) {
	assert(blocksize > 0);
	assert(registers_per_thread > 0);
	assert(shmem_per_block > 0);
	assert(deviceInfo.warpSize > 0 );
	int reg_limit =  deviceInfo.regsPerBlock / (blocksize * registers_per_thread);
	int shm_limit = deviceInfo.sharedMemPerBlock / shmem_per_block ;
	int block_warps = (blocksize+ deviceInfo.warpSize)/deviceInfo.warpSize;
	int total_warps = deviceInfo.maxThreadsPerBlock / deviceInfo.warpSize;
	int warp_limit = block_warps > 0 ? total_warps / block_warps : 0;

	int limit = std::min( warp_limit, std::min( reg_limit , shm_limit ) );

	if(limit == 0)
		$PRINT( "BS: " << blocksize << ", SHM" << shmem_per_block << " -> "
			<< "Limits: reg="	<< reg_limit << ", shm=" << shm_limit 
			<< ", warp=" << warp_limit );

	return limit;
}

bool check_cuda_limits ( int blocksize, int shmem_per_block ){
	return blocks_per_mp(blocksize, shmem_per_block) > 0;
}
