#include <cuda_runtime_api.h>
#include "swarm.h"
#include <vector>
#include <algorithm> // for swap
#include <memory>
#include <iostream>
#include <dlfcn.h>
#include <sstream>
#include <fstream>
#include <valarray>
#include "swarmio.h"

//
// Utilities
//

namespace swarm {
/*!
   \brief  load a configuration file

   @param[out] cfg configuration class
   @param[in] fn file name sting
 */
void load_config(config &cfg, const std::string &fn)
{
	std::ifstream in(fn.c_str());
	if(!in) ERROR("Cannot open configuration file '" + fn + "'.");

	std::string line;
	int iline = 0;
	while(std::getline(in, line))
	{
		iline++;
		line = trim(line);
		if(line.empty()) { continue; }
		if(line[0] == '#') { continue; }

		size_t eqpos = line.find('=');
		if(eqpos == std::string::npos) ERROR("Error on line " + line + ": '=' sign expected.");

		std::string key = trim(line.substr(0, eqpos-1)), val = trim(line.substr(eqpos+1));

		cfg[key] = val;
	}
}

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

/*!
   \brief  Allocate CPU memory for nsys systems of nbod planets each

  @param[in] sys number of systems
  @param[in] nbod number of bodies 
  @param[in] reinitIndices flag for reinitialize indices  
*/
void cpu_ensemble::reset(int nsys, int nbod, bool reinitIndices)	// 
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

//		m_nactive = hostAlloc(m_nactive, 1);
		m_nsys = nsys;
		m_nbod = nbod;
	}

	if(reinitIndices)
	{
		// reinitialize system indices
		for(int i = 0; i != nsys; i++) { m_systemIndices[i] = i; }

		// Set all systems to active
		for(int sys = 0; sys != nsys; sys++) { flags(sys) = 0; }
	}

	m_last_integrator = NULL;
}

/*!
   \brief  Deallocate CPU memory
*/
void cpu_ensemble::free()	
{
//	hostFree(m_nactive); m_nactive = NULL;
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
//	memcpyToHost(m_nactive, src.m_nactive, 1);

	m_last_integrator = src.m_last_integrator;
}

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
//	memcpy(m_nactive, src.m_nactive, 1*sizeof(*m_nactive));

	m_last_integrator = src.m_last_integrator;
}

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
  };
  

void cpu_ensemble::replace_inactive_from(cpu_ensemble &src, const int offset)	// Merge in systems from another cpu_ensemble
{
	int src_sys = 0, dest_sys = 0;
	//	for(int dest_sys=start_sys;dest_sys<nsys();++dest_sys)
	do
	  {
	    //	    do while(is_inactive(dest_sys)&&(dest_sys<nsys())) { ++dest_sys; };
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

// GPU Ensembles ////////////////////////////////
/*!
   \brief GPU Ensemble class 
*/
gpu_ensemble::gpu_ensemble()
{
	construct_base();
	nactive_gpu = NULL;
}

/*!
   \brief GPU Ensemble class 

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
   \brief GPU Ensemble class 

  @param[in] source cpu_ensemble  
*/
gpu_ensemble::gpu_ensemble(const cpu_ensemble &source)
{
	construct_base();
	copy_from(source);
	nactive_gpu = NULL;
}

gpu_ensemble::~gpu_ensemble()
{
	free();
	cudaFree(nactive_gpu);
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

//		cudaMalloc((void**)&m_nactive, sizeof(*m_nactive));
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

		// nactive
//		int tmp3 = nsys;
//		cudaMemcpy(m_nactive, &tmp3, sizeof(tmp3), cudaMemcpyHostToDevice);
	}

	// clear the m_last_integrator field
	m_last_integrator = NULL;
}

/*!
   \brief  Deallocate GPU memory
*/
void gpu_ensemble::free()
{
//	cudaFree(m_nactive); m_nactive = NULL;
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

//int gpu_ensemble::get_nactive() const
//{
//	int nactive;
//	memcpyToHost(&nactive, m_nactive, 1);
//	return nactive;
//}

void debugger_stop()
{
	std::cerr << "Block for debugger here!\n";
}

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
//	memcpyToGPU(m_nactive, src.m_nactive, 1);

	m_last_integrator = src.m_last_integrator;
}

//
// Commonly needed functions
//

/*!
   /brief Calculate totoal energy

   compute the total energy of each system in the ensemble and return it as valarray
  @param[in] ens cpu_ensemble
  @param[out] E total energy
*/
void calc_total_energy(const cpu_ensemble &ens, std::valarray<double> &E)
{
	E.resize(ens.nsys());
	ens.calc_total_energy(&E[0]);
}

/*!
   /brief Ensemble loading support

  @param[in] name string name
  @param[out] ens cpu_ensemble
*/
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

/*!
   /brief Integrator instantiation support

  @param[in] cfg configuration class
*/
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


/*!
   /brief Writer integrator instantiation support

  @param[in] cfg configuration class
*/
writer *writer::create(const std::string &cfg)
{
	std::auto_ptr<writer> w;

	// try loading using a factory function
	void *me = dlopen(NULL, RTLD_LAZY);
	if(me == NULL)
	{
		ERROR(dlerror());
	}

	std::string name;
	std::istringstream ss(cfg);
	if(!(ss >> name))
		ERROR("Empty value for 'output' keyword in config file.");
	std::string factory_name = "create_writer_" + name;

	writerFactory_t factory = (writerFactory_t)dlsym(me, factory_name.c_str());
	if(factory)
	{
		std::string wcfg;
		getline(ss, wcfg);
		wcfg = trim(wcfg);
		w.reset(factory(wcfg));
	}
	else
	{
		ERROR("Writer " + name + " unknown (" + dlerror() + ").");
	}

	return w.release();
}

/*!
   /brief Find best factorization 

   Find the dimensions (bx,by) of a 2D grid of blocks that has as close to nblocks blocks as possible
  @param[out] bx
  @param[out] by
  @param[in] nblocks
*/
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

/*!
   /brief Configur grid

   Given a total number of threads, their memory requirements, and the
   number of threadsPerBlock, compute the optimal allowable grid dimensions.
   Returns false if the requested number of threads are impossible to fit to
   shared memory.

  @param[out] gridDim
  @param[in] threadsPerBlock
  @param[in] nthreads
  @param[in] dynShmemPerThread
  @param[in] staticShmemPerBlcok
 */
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


} // end namespace swarm
