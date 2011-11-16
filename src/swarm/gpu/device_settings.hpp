


int shmem_per_mpu()  {
	  int cuda_dev_id = -1;
	  cuxErrCheck( cudaGetDevice(&cuda_dev_id) );
	  assert(cuda_dev_id>=0);
	  cudaDeviceProp deviceProp;
	  cuxErrCheck( cudaGetDeviceProperties(&deviceProp, cuda_dev_id) );
	  return deviceProp.sharedMemPerBlock;
}


/**
 *  This function should optimize based on 3 properties of the device and 1 compile time value
 *   - Shared Memory per MPU
 *   - Maximum block size
 *   - Device Warp size
 *   - Ensemble chunk size
 */
int optimized_system_per_block(const int thread_per_system,const int shmem_per_system)  {
	//	  const int nbod = _hens.nbod();
	const int max_shmem = shmem_per_mpu();
	//	  return 16;
	// WARNING: Assumes that integrator does not use shared memory beyond what is used by standard optimized Gravitation class.  
	// Can integrators like hermite_adap override this?
	const int max_system_per_block_with_two_blocks = max_shmem/(2*shmem_per_system);
	const int default_system_per_block = std::max(1,max_system_per_block_with_two_blocks);
	// WARNING: Does not account for larger memory usage due to coalesced arrys.  Is this the reason I had to set CHUNK_SIZE=1 in gravitation?
	std::cerr << "# default_system_per_block = " << default_system_per_block << " nbod = " << nbod << " max_shmem= " << max_shmem << " shmem_per_system(nbod)= " << shmem_per_system(nbod) << "\n";
	return default_system_per_block;
}

