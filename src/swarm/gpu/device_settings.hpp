
void select_cuda_device(int dev);
void print_device_information();
int blocks_per_mp( int blocksize, int shmem_per_block ) ;
bool check_cuda_limits ( int blocksize, int shmem_per_block );
int optimized_system_per_block(int chunk_size, int thread_per_system
		, int shmem_per_system);
void set_more_cache();
