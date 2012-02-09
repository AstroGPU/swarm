#include "swarm/common.hpp"
#include "device_settings.hpp"

// Setting prefered cache size the following function
// can do this but the problem is that since our kernel is really
// a template functin we don't know the exact name for it
set_prefere_shared_memory(const char* function_name){
	cudaFuncSetCacheConfig( name, cudaFuncCachePreferShared);
}

const int registers_per_thread = 64;
cudaDeviceProp deviceInfo;

void select_cuda_device(int dev) {
	int devcnt; cudaErrCheck( cudaGetDeviceCount(&devcnt) );
	if( dev >= 0 && dev < devcnt )
		cudaErrCheck( cudaSetDevice(dev) );
	else
		std::cerr << "Cannot select the CUDA device. GPU integrators are disabled" << std::endl;

	  cudaErrCheck( cudaGetDeviceProperties(&deviceInfo, dev) );
}



void print_device_information(){
	  std::cout << "Device:\t"  << deviceInfo.name  
		  << " Compute Capabality: " << deviceInfo.computeMode <<   "\n"
		  << "Global Memory:\t" << deviceInfo.totalGlobalMem/double(1<<30)     << "GB\n"
		  << "Shared Memory\t"  << deviceInfo.sharedMemPerBlock/1024  << "KB\n"
		  << "Max Blocksize\t"  << deviceInfo.maxThreadsPerBlock << "\n"
		  << "Warp Size    \t"  << deviceInfo.warpSize   << "\n"
		  << "Registers/MP \t"  << deviceInfo.regsPerBlock       << "\n"
		  ;
	  
}


int blocks_per_mp( int blocksize, int shmem_per_block ) {
	int reg_limit =  deviceInfo.regsPerBlock / (blocksize * registers_per_thread);
	int shm_limit = deviceInfo.sharedMemPerBlock / shmem_per_block ;
	int block_warps = (blocksize+ deviceInfo.warpSize)/deviceInfo.warpSize;
	int total_warps = deviceInfo.maxThreadsPerBlock / deviceInfo.warpSize;
	int warp_limit = total_warps / block_warps;

	int limit = std::min( warp_limit, std::min( reg_limit , shm_limit ) );

	if(limit == 0)
		$PRINT( reg_limit << ", " << shm_limit << ", " << warp_limit );

	return limit;
}

bool check_cuda_limits ( int blocksize, int shmem_per_block ){
	return blocks_per_mp(blocksize, shmem_per_block) > 0;
}
