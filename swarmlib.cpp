#include <cuda_runtime_api.h>
#include "swarm.h"
#include "integrators.h"
#include <vector>
#include <iostream>

//
// Utilities
//

void die(const std::string &msg)
{
	std::cerr << msg << "\n";
	abort();
}

template<typename T>
inline void memcpyToGPU(T *dest, const T *src, int nelem)
{
	cudaMemcpy(dest, src, nelem*sizeof(T), cudaMemcpyHostToDevice);
}

template<typename T>
inline void memcpyToHost(T *dest, const T *src, int nelem)
{
	cudaMemcpy(dest, src, nelem*sizeof(T), cudaMemcpyDeviceToHost);
}

//
// Ensemble class plumbing (mostly memory management)
//

// CPU Ensembles ////////////////////////////////
cpu_ensemble::cpu_ensemble()
{
	construct_base();
}

cpu_ensemble::cpu_ensemble(int nsys, int nbod)
{
	construct_base();
	reset(nsys, nbod);
}

cpu_ensemble::cpu_ensemble(const gpu_ensemble &source)
{
	construct_base();
	copy_from(source);
}

void cpu_ensemble::reset(int nsys, int nbod, bool reinitIndices)	// Allocate CPU memory for nsys systems of nbod planets each
{
	// do we need to realloc?
	if(m_nsys != nsys || m_nbod != nbod)
	{
		m_T = (float*)realloc(m_T, nsys*sizeof(*m_T));
		m_x = (double*)realloc(m_x, nsys*nbod*sizeof(*m_x));
		m_y = (double*)realloc(m_y, nsys*nbod*sizeof(*m_y));
		m_z = (double*)realloc(m_z, nsys*nbod*sizeof(*m_z));
		m_vx = (double*)realloc(m_vx, nsys*nbod*sizeof(*m_vx));
		m_vy = (double*)realloc(m_vy, nsys*nbod*sizeof(*m_vy));
		m_vz = (double*)realloc(m_vz, nsys*nbod*sizeof(*m_vz));
		m_m = (float*)realloc(m_m, nsys*nbod*sizeof(*m_m));
		m_active = (int*)realloc(m_active, nsys*sizeof(*m_active));
		m_systemIndices = (int*)realloc(m_systemIndices, nsys*sizeof(*m_systemIndices));

		m_nsys = nsys;
		m_nbod = nbod;
	}

	if(reinitIndices)
	{
		// reinitialize system indices
		for(int i = 0; i != nsys; i++) { m_systemIndices[i] = i; }
	}

	m_last_integrator = NULL;
}

void cpu_ensemble::free()			// Deallocate CPU memory
{
	::free(m_T); m_T = NULL;
	::free(m_x); m_x = NULL;
	::free(m_y); m_y = NULL;
	::free(m_z); m_z = NULL;
	::free(m_vx); m_vx = NULL;
	::free(m_vy); m_vy = NULL;
	::free(m_vz); m_vz = NULL;
	::free(m_m); m_m = NULL;
	::free(m_active); m_active = NULL;
	::free(m_systemIndices); m_systemIndices = NULL;
}

void cpu_ensemble::copy_from(const gpu_ensemble &src)	// Copy the data from the GPU
{
	reset(src.nsys(), src.nbod(), false);

	// low-level copy from host to device memory
	memcpyToHost(m_T, src.m_T, m_nsys);
	memcpyToHost(m_x, src.m_x, m_nbod*m_nsys);
	memcpyToHost(m_y, src.m_y, m_nbod*m_nsys);
	memcpyToHost(m_z, src.m_z, m_nbod*m_nsys);
	memcpyToHost(m_vx, src.m_vx, m_nbod*m_nsys);
	memcpyToHost(m_vy, src.m_vy, m_nbod*m_nsys);
	memcpyToHost(m_vz, src.m_vz, m_nbod*m_nsys);
	memcpyToHost(m_m, src.m_m, m_nsys);
	memcpyToHost(m_active, src.m_active, m_nsys);
	memcpyToHost(m_systemIndices, src.m_systemIndices, m_nsys);

	m_last_integrator = src.m_last_integrator;
	m_nactive = src.m_nactive;
}

// GPU Ensembles ////////////////////////////////
gpu_ensemble::gpu_ensemble()
{
	construct_base();
}

gpu_ensemble::gpu_ensemble(int nsys, int nbod)
{
	construct_base();
	reset(nsys, nbod);
}

gpu_ensemble::gpu_ensemble(const cpu_ensemble &source)
{
	construct_base();
	copy_from(source);
}

void gpu_ensemble::reset(int nsys, int nbod, bool reinitIndices)	// Allocate CPU memory for nsys systems of nbod planets each
{
	// do we need to realloc?
	if(m_nsys != nsys || m_nbod != nbod)
	{
		free();

		cudaMalloc((void**)&m_T, nsys*sizeof(*m_T));
		cudaMalloc((void**)&m_x, nsys*nbod*sizeof(*m_x));
		cudaMalloc((void**)&m_y, nsys*nbod*sizeof(*m_y));
		cudaMalloc((void**)&m_z, nsys*nbod*sizeof(*m_z));
		cudaMalloc((void**)&m_vx, nsys*nbod*sizeof(*m_vx));
		cudaMalloc((void**)&m_vy, nsys*nbod*sizeof(*m_vy));
		cudaMalloc((void**)&m_vz, nsys*nbod*sizeof(*m_vz));
		cudaMalloc((void**)&m_m, nsys*nbod*sizeof(*m_m));
		cudaMalloc((void**)&m_active, nsys*sizeof(*m_active));
		cudaMalloc((void**)&m_systemIndices, nsys*sizeof(*m_systemIndices));

		m_nsys = nsys;
		m_nbod = nbod;
	}

	if(reinitIndices)
	{
		// reinitialize system indices
		std::vector<int> tmp(nsys);
		for(int i = 0; i != nsys; i++) { tmp[i] = i; }
		cudaMemcpy(m_systemIndices, &tmp[0], nsys, cudaMemcpyHostToDevice);
	}

	// clear the m_last_integrator field
	m_last_integrator = NULL;
}

void gpu_ensemble::free()			// Deallocate CPU memory
{
	cudaFree(m_T); m_T = NULL;
	cudaFree(m_x); m_x = NULL;
	cudaFree(m_y); m_y = NULL;
	cudaFree(m_z); m_z = NULL;
	cudaFree(m_vx); m_vx = NULL;
	cudaFree(m_vy); m_vy = NULL;
	cudaFree(m_vz); m_vz = NULL;
	cudaFree(m_m); m_m = NULL;
	cudaFree(m_active); m_active = NULL;
	cudaFree(m_systemIndices); m_systemIndices = NULL;
}

void gpu_ensemble::copy_from(const cpu_ensemble &src)	// Copy the data from the GPU
{
	reset(src.nsys(), src.nbod(), false);
	
	// low-level copy from host to device memory
	memcpyToGPU(m_T, src.m_T, m_nsys);
	memcpyToGPU(m_x, src.m_x, m_nbod*m_nsys);
	memcpyToGPU(m_y, src.m_y, m_nbod*m_nsys);
	memcpyToGPU(m_z, src.m_z, m_nbod*m_nsys);
	memcpyToGPU(m_vx, src.m_vx, m_nbod*m_nsys);
	memcpyToGPU(m_vy, src.m_vy, m_nbod*m_nsys);
	memcpyToGPU(m_vz, src.m_vz, m_nbod*m_nsys);
	memcpyToGPU(m_m, src.m_m, m_nsys);
	memcpyToGPU(m_active, src.m_active, m_nsys);
	memcpyToGPU(m_systemIndices, src.m_systemIndices, m_nsys);

	m_last_integrator = src.m_last_integrator;
	m_nactive = src.m_nactive;
}


//
// Integrator support
//

integrator *create_integrator(const std::string &type, const config &cfg)
{
	if(type == "gpu_euler") { return new gpu_euler_integrator(cfg); }
	else return NULL;
}

gpu_euler_integrator::gpu_euler_integrator(const config &cfg)
{
	assert(cfg.count("h"));
	h = atof(cfg.at("h").c_str());
}


//
// Find the dimensions (bx,by) of a 2D grid of blocks that 
// has as close to nblocks blocks as possible
//
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

//
// Given a total number of threads, their memory requirements, and the
// number of threadsPerBlock, compute the optimal allowable grid dimensions.
// Returns false if the requested number of threads are impossible to fit to
// shared memory.
//
bool configure_grid(dim3 &gridDim, int &threadsPerBlock, int nthreads, int dynShmemPerThread, int staticShmemPerBlock)
{
	const int shmemPerMP =  16384;
	threadsPerBlock = 192; // HACK: should compute this dynamically, based on memory requirements

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

	std::cerr << "Grid: Threads requested = " << nthreads << " with shmem/thread = " << dynShmemPerThread << " and shmem/blk = " << staticShmemPerBlock << "\n";
	std::cerr << "      Grid configured as (" << gridDim.x << ", " << gridDim.y << ", " << gridDim.z <<") array of blocks with " << threadsPerBlock << " threads per block.\n";
	std::cerr << "      Total threads to execute = " << nthreadsEx << "\n";

	return true;
}
