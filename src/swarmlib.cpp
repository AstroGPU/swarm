#include <cuda_runtime_api.h>
#include "swarm.h"
#include <vector>
#include <iostream>
#include <dlfcn.h>
#include <sstream>
#include <fstream>
#include <valarray>
#include "swarmio.h"

//
// Utilities
//

void die(const std::string &msg)
{
	std::cerr << msg << "\n";
	abort();
}

void trim(std::string& str)
{
	std::string::size_type pos = str.find_last_not_of(" \t");
	if (pos != std::string::npos)
	{
		str.erase(pos + 1);
		pos = str.find_first_not_of(" \t");
		if (pos != std::string::npos) str.erase(0, pos);
	}
	else str.erase(str.begin(), str.end());
}

// load a configuration file
void load_config(config &cfg, const std::string &fn)
{
	std::ifstream in(fn.c_str());
	if(!in) ERROR("Cannot open configuration file '" + fn + "'.");

	std::string line;
	int iline = 0;
	while(std::getline(in, line))
	{
		iline++;
		trim(line);
		if(line.empty()) { continue; }
		if(line[0] == '#') { continue; }

		size_t eqpos = line.find(' ');
		if(eqpos == std::string::npos) ERROR("Error on line " + str(line) + ": '=' sign expected.");

		std::string key = line.substr(0, eqpos), val = line.substr(eqpos+2);
		trim(key); trim(val);

		cfg[key] = val;
	}
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
		m_T = hostAlloc(m_T, nsys);
		m_xyz = hostAlloc(m_xyz, 3*nsys*nbod);
		m_vxyz = hostAlloc(m_vxyz, 3*nsys*nbod);
		m_m = hostAlloc(m_m, nsys*nbod);
		m_active = hostAlloc(m_active, nsys);
		m_systemIndices = hostAlloc(m_systemIndices, nsys);

		m_nsys = nsys;
		m_nbod = nbod;
	}

	if(reinitIndices)
	{
		// reinitialize system indices
		for(int i = 0; i != nsys; i++) { m_systemIndices[i] = i; }

		// Set all systems to active
		for(int sys = 0; sys != nsys; sys++) { active(sys) = true; }
	}

	m_last_integrator = NULL;
}

void cpu_ensemble::free()			// Deallocate CPU memory
{
	hostFree(m_T); m_T = NULL;
	hostFree(m_xyz); m_xyz = NULL;
	hostFree(m_vxyz); m_vxyz = NULL;
	hostFree(m_m); m_m = NULL;
	hostFree(m_active); m_active = NULL;
	hostFree(m_systemIndices); m_systemIndices = NULL;
}

void cpu_ensemble::copy_from(const gpu_ensemble &src)	// Copy the data from the GPU
{
	reset(src.nsys(), src.nbod(), false);

	// low-level copy from host to device memory
	memcpyToHost(m_T, src.m_T, m_nsys);
	memcpyToHost(m_xyz, src.m_xyz, 3*m_nbod*m_nsys);
	memcpyToHost(m_vxyz, src.m_vxyz, 3*m_nbod*m_nsys);
	memcpyToHost(m_m, src.m_m, m_nbod*m_nsys);
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
		cudaMalloc((void**)&m_xyz, 3*nsys*nbod*sizeof(*m_xyz));
		cudaMalloc((void**)&m_vxyz, 3*nsys*nbod*sizeof(*m_vxyz));
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
		cudaMemcpy(m_systemIndices, &tmp[0], tmp.size()*sizeof(tmp[0]), cudaMemcpyHostToDevice);

		// Set all systems active
		std::valarray<int> tmp2(nsys);
		tmp2 = 1;
		cudaMemcpy(m_active, &tmp2[0], tmp2.size()*sizeof(tmp2[0]), cudaMemcpyHostToDevice);
	}

	// clear the m_last_integrator field
	m_last_integrator = NULL;
}

void gpu_ensemble::free()			// Deallocate CPU memory
{
	cudaFree(m_T); m_T = NULL;
	cudaFree(m_xyz); m_xyz = NULL;
	cudaFree(m_vxyz); m_vxyz = NULL;
	cudaFree(m_m); m_m = NULL;
	cudaFree(m_active); m_active = NULL;
	cudaFree(m_systemIndices); m_systemIndices = NULL;
}

void gpu_ensemble::copy_from(const cpu_ensemble &src)	// Copy the data from the GPU
{
	reset(src.nsys(), src.nbod(), false);
	
	// low-level copy from host to device memory
	memcpyToGPU(m_T, src.m_T, m_nsys);
	memcpyToGPU(m_xyz, src.m_xyz, 3*m_nbod*m_nsys);
	memcpyToGPU(m_vxyz, src.m_vxyz, 3*m_nbod*m_nsys);
	memcpyToGPU(m_m, src.m_m, m_nbod*m_nsys);
	memcpyToGPU(m_active, src.m_active, m_nsys);
	memcpyToGPU(m_systemIndices, src.m_systemIndices, m_nsys);

	m_last_integrator = src.m_last_integrator;
	m_nactive = src.m_nactive;
}

//
// Commonly needed functions
//

// compute the total energy of each system in the ensemble and return it as valarray
void calc_total_energy(const cpu_ensemble &ens, std::valarray<double> &E)
{
	E.resize(ens.nsys());
	ens.calc_total_energy(&E[0]);
}

//
// Ensemble loading support
//
void load_ensemble(const std::string &name, cpu_ensemble &ens)
{
	// Assume filenames are of the form "$name.xx" where xx is a sequence of
	// numbers 0..(NSYS-1).
	//
	// Initial conditions file format:
	// <nbod>
	// <m1> <x1> <y1> <z1> <vx1> <vy1> <vz1>
	// <m2> <x2> <y2> <z2> <vx2> <vy2> <vz2>
	// ... (through nbod) ..
	//

	// determine the number of systems
	int nsys;
	for(nsys = 0; true; nsys++)
	{
		std::ostringstream ss; ss << nsys;
		std::string fn = name + "." + ss.str();
		if(access(fn.c_str(), R_OK) != 0) { break; }
	}
	if(nsys == 0) ERROR("Failed to find the datafiles for ensemble " + name + " (cannot access " + name + ".0)");

	int nbod = -1;
	for(int i = 0; i != nsys; i++)
	{
		std::ostringstream ss; ss << i;
		std::string fn = name + "." + ss.str();
		std::ifstream in(fn.c_str());

		// load the number of bodies
		int nbod1;
		in >> nbod1;
		if(nbod == -1)
		{
			nbod = nbod1;
			ens.reset(nsys, nbod);
		}
		else if(nbod != nbod1)
		{
			std::ostringstream err;
			err << "The number of bodies must be the same for all systems. Found nbod=" << nbod1 << " in " << fn
			    << " while expecting nbod=" << nbod << ".";
			ERROR(err.str());
		}

		// set the initial time to 0 (NOTE: should we also load this from the config file?)
		ens.time(i) = 0.;

		// load the planets
		double m, x, y, z, vx, vy, vz;
		for(int j = 0; j != nbod; j++)
		{
			if(!(in >> m >> x >> y >> z >> vx >> vy >> vz))
			{
				ERROR("Error loading bodies from " + fn + ".");
			}
			ens.set_body(i, j, m, x, y, z, vx, vy, vz);
		}
	}

	std::cerr << "Loaded " << nsys << " systems of " << nbod << " bodies each.\n";
}

//
// Integrator instantiation support
//

integrator *integrator::create(const config &cfg)
{
	std::auto_ptr<integrator> integ;

	// try loading using a factory function
	void *me = dlopen(NULL, RTLD_LAZY);
	if(me == NULL)
	{
		ERROR(dlerror());
	}

	if(!cfg.count("integrator")) ERROR("Integrator type has to be chosen with 'integrator=xxx' keyword");
	
	std::string name = cfg.at("integrator");
	std::string factory_name = "create_" + name;

	integratorFactory_t factory = (integratorFactory_t)dlsym(me, factory_name.c_str());
	if(factory)
	{
		integ.reset(factory(cfg));
	}
	else
	{
		ERROR("Integrator " + name + " unknown (" + dlerror() + ").");
	}

	return integ.release();
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

#if 0
	std::cerr << "+ Grid configuration =========================\n";
	std::cerr << "      Threads requested = " << nthreads << " with shmem/thread = " << dynShmemPerThread << " and shmem/blk = " << staticShmemPerBlock << "\n";
	std::cerr << "      Grid configured as (" << gridDim.x << ", " << gridDim.y << ", " << gridDim.z <<") array of blocks with " << threadsPerBlock << " threads per block.\n";
	std::cerr << "      Total threads to execute = " << nthreadsEx << "\n";
	std::cerr << "- Grid configuration =========================\n";
#else
	std::cerr << "Kernel exec. config: (" << gridDim.x << ", " << gridDim.y << ", " << gridDim.z <<") x " << threadsPerBlock << " thr/blk (" << nthreadsEx << " thr total; " << nthreads << " thr needed)\n";
#endif
	return true;
}

//
// I/O and snapshotting functions
//
ens_writer::ens_writer(const std::string &fn_)
	: fn(fn_), out(fn.c_str()), bout(out)
{
	if(!out) ERROR("Problem opening output file '" + fn + "'");
}

ens_writer &ens_writer::operator <<(const cpu_ensemble &ens)
{
	bout << ens.nsys() << ens.nbod() << ens.nactive();

	for(int sysID=0; sysID != ens.nsys(); sysID++)
	{
		int sys = ens.index_of_system(sysID);

		bout << ens.time(sys);
		bout << ens.active(sys);

		for(int bod = 0; bod != ens.nbod(); bod++)
		{
			bout << ens.mass(sys, bod);
			bout << ens.x(sys, bod) << ens.y(sys, bod) << ens.z(sys, bod);
			bout << ens.vx(sys, bod) << ens.vy(sys, bod) << ens.vz(sys, bod);
		}
	}

	return *this;
}

ens_reader::ens_reader(const std::string &fn_)
	: fn(fn_), in(fn.c_str()), bin(in)
{
	if(!in) ERROR("Problem opening input file '" + fn + "'");
}

ens_reader &ens_reader::operator >>(cpu_ensemble &ens)
{
	int nsys, nbod, nactive;
	if(!(bin >> nsys))
	{
		//if(in.eof())
			return *this;
		//ERROR("Data file " + fn + " corrupted.");
	}

	if(!(bin >> nbod >> nactive))
		ERROR("Data file " + fn + " corrupted.");

	ens.reset(nsys, nbod);

	for(int sys=0; sys != ens.nsys(); sys++)
	{
		bin >> ens.time(sys);
		bin >> ens.active(sys);

		for(int bod = 0; bod != ens.nbod(); bod++)
		{
			bin >> ens.mass(sys, bod);
			bin >> ens.x(sys, bod) >> ens.y(sys, bod) >> ens.z(sys, bod);
			bin >> ens.vx(sys, bod) >> ens.vy(sys, bod) >> ens.vz(sys, bod);
		}

		if(!bin) ERROR("Data file " + fn + " corrupted.");
	}

	return *this;
}
