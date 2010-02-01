#ifndef __swarm_h
#define __swarm_h

#include <stdexcept>
#include <string>
#include <cstring>
#include <map>
#include <cassert>
#include <cmath>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cux/cux.h>

#include "stopwatch.h"

class integrator;
class ieventstream;
class writer;
class gpu_ensemble;
class cpu_ensemble;

/**
	\brief Unrecoverable error exception.

	Throw an instance of this class to indicate an unrecoverable error
	was encountered. Do not throw it directly, but through the use of ERROR() macro.
*/
class swarm_error : public std::runtime_error
{
public:
	swarm_error(const std::string &msg) : std::runtime_error(msg) {}
	virtual ~swarm_error() throw() {};
};

#ifndef THROW_IS_ABORT
	#define ERROR(msg) throw swarm_error(msg);
#else
	#include <cassert>
	#include <cstring>
	#define ERROR(msg) { fprintf(stderr, "%s\n", std::string(msg).c_str()); abort(); }
#endif

/**
	\brief Ensemble data storage base class

	Ensemble base class. Defines the minimum subset of data commonly needed
	by all integrators (as well as the user). Must not be used directly by
	the user; use cpu_ensemble and gpu_ensemble instead.

	Note that it has no constructors/destructors/virtuals, to allow its
	instantiation in __constant__ memory on the GPU.
*/
typedef double real_time;
typedef float  real_mass;
typedef double real_pos;
typedef double real_vel;

class ensemble
{
	public:
		enum { INACTIVE = 0x01 };

	protected:
		// number of active (currently integrating) systems, total number of systems
		// and number of bodies per system
		int m_nsys, m_nbod;

		// m_nsys wide array
		real_time *m_T;
		real_time *m_Tend;

		// m_nsys*2 wide arrays
		real_time *m_Toutput;

		// m_nsys*m_nbod*3 wide arrays
		real_pos  *m_xyz;
		real_vel  *m_vxyz;

		// m_nsys*m_nbod wide arrays
		real_mass *m_m;

		// m_nsys wide arrays
		int	*m_flags;			// flags about a given system. Bit 0 is the inactivity flag (m_flags & 0x01 ? inactive : acctive). Others can be user-defined.
		int	*m_systemIndices;		// map from original systemID to sys index in m_* arrays
		integrator *m_last_integrator;

		// scalars
//		int	*m_nactive;			// number of active systems, sum of !(m_flags & 0x01). Computed on GPU on exit from kernel.

	public:
		// DEPRECATED: For Young In's code ONLY!!!
		double *xyz()  { return m_xyz; }
		double *vxyz() { return m_vxyz; }
		const double *xyz()  const { return m_xyz; }
		const double *vxyz() const { return m_vxyz; }

	protected:
		void construct_base()
		{
			memset(this, 0, sizeof(*this));
		}

	public:
		// non-const versions
		__host__ __device__ double&   time(int sys) { return m_T[sys]; }
		__host__ __device__ double&   time_end(int sys) { return m_Tend[sys]; }
		__host__ __device__ double&   time_output(int sys, int k) { return m_Toutput[k*m_nsys + sys]; }
	
		__host__ __device__ double&  x(int sys, int bod) { return m_xyz[bod*m_nsys + sys]; }
		__host__ __device__ double&  y(int sys, int bod) { return m_xyz[m_nbod*m_nsys + bod*m_nsys + sys]; }
		__host__ __device__ double&  z(int sys, int bod) { return m_xyz[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

		__host__ __device__ double& vx(int sys, int bod) { return m_vxyz[bod*m_nsys + sys]; }
		__host__ __device__ double& vy(int sys, int bod) { return m_vxyz[m_nbod*m_nsys + bod*m_nsys + sys]; }
		__host__ __device__ double& vz(int sys, int bod) { return m_vxyz[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

		__host__ __device__ float& mass(int sys, int bod)   { return m_m[bod*m_nsys + sys]; }

		__host__ __device__ int& flags(int sys)	{ return m_flags[sys]; }

//		__host__ __device__ int& nactive() { return *m_nactive; }
		__host__ __device__ int& nsys() { return m_nsys; }
		__host__ __device__ int& nbod() { return m_nbod; }

		__host__ __device__ int& index_of_system(int sysId) { return m_systemIndices[sysId]; }


		// const versions
		__host__ __device__ double time(int sys) const { return m_T[sys]; }
		__host__ __device__ double time_end(int sys) const { return m_Tend[sys]; }
		__host__ __device__ double time_output(int sys, int k) const { return m_Toutput[k*m_nsys + sys]; }
	
		__host__ __device__ double  x(int sys, int bod) const { return m_xyz[bod*m_nsys + sys]; }
		__host__ __device__ double  y(int sys, int bod) const { return m_xyz[m_nbod*m_nsys + bod*m_nsys + sys]; }
		__host__ __device__ double  z(int sys, int bod) const { return m_xyz[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

		__host__ __device__ double vx(int sys, int bod) const { return m_vxyz[bod*m_nsys + sys]; }
		__host__ __device__ double vy(int sys, int bod) const { return m_vxyz[m_nbod*m_nsys + bod*m_nsys + sys]; }
		__host__ __device__ double vz(int sys, int bod) const { return m_vxyz[m_nbod*m_nsys*2 + bod*m_nsys + sys]; }

		__host__ __device__ float mass(int sys, int bod) const { return m_m[bod*m_nsys + sys]; }

		__host__ __device__ int flags(int sys)		const { return m_flags[sys]; }

//		__host__ __device__ int nactive() const { return *m_nactive; }
		__host__ __device__ int nsys() const { return m_nsys; }
		__host__ __device__ int nbod() const { return m_nbod; }

		__host__ __device__ int index_of_system(int sysId) const { return m_systemIndices[sysId]; }



		// convenience
		__host__ __device__ int active(int sys)		const { return m_flags[sys] ^ ~ensemble::INACTIVE; }

		__host__ __device__ void set_body(int sys, int bod,  float m, double x, double y, double z, double vx, double vy, double vz)
		{
			int idx = bod*m_nsys + sys;

			m_m[idx]   =  m;
			m_xyz[idx] =  x; m_xyz[m_nbod*m_nsys + idx] =  y;  m_xyz[m_nbod*m_nsys*2 + idx] =  z;
			m_vxyz[idx]  = vx;  m_vxyz[m_nbod*m_nsys + idx] = vy;   m_vxyz[m_nbod*m_nsys*2 + idx] = vz;
		}

		__host__ __device__ void get_body(int sys, int bod, float &m, double &x, double &y, double &z, double &vx, double &vy, double &vz) const
		{
			int idx = bod*m_nsys + sys;
			
			 m = m_m[idx];
			 x = m_xyz[idx];   y = m_xyz[m_nbod*m_nsys + idx];   z = m_xyz[m_nbod*m_nsys*2 + idx];
			vx = m_vxyz[idx]; vy = m_vxyz[m_nbod*m_nsys + idx]; vz = m_vxyz[m_nbod*m_nsys*2 + idx];
		}

		// utilities
		__host__ __device__ double calc_total_energy(int sys) const
		{
			double E = 0.;
			for (int bod1 = 0; bod1 != nbod(); bod1++)
			{
				float m1; double x1[3], v1[3];
				get_body(sys, bod1, m1, x1[0], x1[1], x1[2], v1[0], v1[1], v1[2]);
				E += 0.5 * m1 * (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);

				for (int bod2 = 0; bod2 < bod1; bod2++)
				{
					float m2; double x2[3], v2[3];
					get_body(sys, bod2, m2, x2[0], x2[1], x2[2], v2[0], v2[1], v2[2]);
					double dist = sqrt((x2[0] - x1[0]) * (x2[0] - x1[0]) + (x2[1] - x1[1]) * (x2[1] - x1[1]) + (x2[2] - x1[2]) * (x2[2] - x1[2]));

					E -= m1 * m2 / dist;
				}
			}
			return E;
		}

		__host__ __device__ void calc_total_energy(double *E) const
		{
			for (int sys = 0; sys != nsys(); sys++)
			{
				E[sys] = calc_total_energy(sys);
			}
		}

	public:
		integrator *last_integrator() { return m_last_integrator; }
		void set_last_integrator(integrator *li) { m_last_integrator = li; }

	private:
		friend class cpu_ensemble;
		friend class gpu_ensemble;
};

/**
	\brief Specialization of ensemble, with data stored on the host (CPU)

	Use this class to create/read/write ensembles in host memory.
*/
class cpu_ensemble : public ensemble
{
	public:
		void reset(int nsys, int nbod, bool reinitIndices = true);	// Allocate CPU memory for nsys systems of nbod planets each
		void free();						// Deallocate CPU memory

	public:
		void copy_from(const gpu_ensemble &source);	// Copy the data from the GPU
		
	public:
		cpu_ensemble();					// instantiate an empty ensemble
		cpu_ensemble(int nsys, int nbod);		// instantiate an ensemble with room for nsys systems with nbod each
		cpu_ensemble(const gpu_ensemble &source);	// instantiate a copy of source ensemble

		~cpu_ensemble() { free(); }

	private:
		cpu_ensemble &operator=(const cpu_ensemble&);	// disallow copying
};

/**
	\brief Specialization of ensemble, with data stored on the device (GPU)

	Use this class to create/read/write ensembles in GPU memory. Typically, you
	will want to create a cpu_ensemble, populate it on the host, use gpu_ensemble's
	copy constructor to upload the data to GPU, and pass the gpu_ensemble structure
	to an integrator that will advance the ensemble.
*/
class gpu_ensemble : public ensemble
{
	protected:
		int *nactive_gpu;	// temp variable for get_nactive()

	public:
		void reset(int nsys, int nbod, bool reinitIndices = true);	// Allocate GPU memory for nsys systems of nbod planets each
		void free();							// Deallocate CPU memory

	public:
		void copy_from(const cpu_ensemble &source);	// Copy the data from the CPU
		int get_nactive() const;			// Download and return the number of active systems

	public:
		gpu_ensemble();					// instantiate an empty ensemble
		gpu_ensemble(int nsys, int nbod);		// instantiate an ensemble with room for nsys systems with nbod each
		gpu_ensemble(const cpu_ensemble &source);	// instantiate a copy of the source ensemble

		~gpu_ensemble();
};

typedef std::map<std::string, std::string> config;

/// Load ensemble residing in files "name.XXX" where XXX \elem [0,nsys)
void load_ensemble(const std::string &name, cpu_ensemble &ens);
/// Load configuration from file fn
void load_config(config &cfg, const std::string &fn);

#ifndef __CUDACC__ // CUDA 2.2 C++ bug workaround
#include <sstream>

// get a configuration value for 'key', throwing an error if it doesn't exist
// NOTE: heavy (unoptimized) function, use sparingly
template<typename T>
void get_config(T &val, const config &cfg, const std::string &key)
{
	if(!cfg.count(key)) { ERROR("Configuration key '" + key + "' missing."); }
	std::istringstream ss(cfg.at(key));
	ss >> val;
}
#endif

template<typename T>
inline void memcpyToGPU(T *dest, const T *src, int nelem = 1)
{
	cuxErrCheck( cudaMemcpy(dest, src, nelem*sizeof(T), cudaMemcpyHostToDevice) );
}

template<typename T>
inline void memcpyToHost(T *dest, const T *src, int nelem = 1)
{
	cuxErrCheck( cudaMemcpy(dest, src, nelem*sizeof(T), cudaMemcpyDeviceToHost) );
}

/**
	\brief Abstract output writer interface

	The method process() is called whenever the GPU output buffers are filled,
	and should proces/store the accumulated output data and logs (usually by
	writing them out to a file).
*/
class writer
{
	public:
		virtual void process(ieventstream &es) = 0;
		virtual ~writer() {};	// has to be here to ensure the derived class' destructor is called (if it exists)

	public:
		// Integrator factory functions (and supporting typedefs)
		static writer *create(const std::string &cfg);

	protected:
		writer() {};		// hide the constructor.and force integrator instantiation with integrator::create
};
typedef writer *(*writerFactory_t)(const std::string &cfg);

/**
	\brief Abstract integrator interface

	Specific integrators must derive from this class (as well as define a
	create_XXXX(const config &cfg) factory function (where XXXX is the name
	of the integrator)). Instantiate them via a call to integrator::create(),
	where cfg["integrator"] is expected to contain the name of the integrator
	you wish to instantiate.
*/
class integrator
{
	public:
		virtual void integrate(gpu_ensemble &ens, double T)	// for GPU based integrators
			{ ERROR("Execution on GPU not supported by this implementation"); }
		virtual void integrate(cpu_ensemble &ens, double T)	// for CPU based integrators
			{ ERROR("Execution on GPU not supported by this implementation"); }

		virtual ~integrator() {};	// has to be here to ensure the derived class' destructor is called (if it exists)

	public:
		// Integrator factory functions (and supporting typedefs)
		static integrator *create(const config &cfg);

	protected:
		integrator() {};		// hide the constructor.and force integrator instantiation with integrator::create
};

typedef integrator *(*integratorFactory_t)(const config &cfg);

// configure grid for nthreads intependent threads, each requiring dynShmemPerThread of shared memory, and
// with each block needing staticShmemPerBlock of shared memory (usually to pass kernel invocation args.)
// The default for staticShmemPerBlock is reasonable (for small kernels), but not necessarily optimal.
bool configure_grid(dim3 &gridDim, int threadsPerBlock, int nthreads, int dynShmemPerThread = 0, int staticShmemPerBlock = 128);

// Typesafe de-allocator (convenience)
template<typename T>
void hostFree(T* var, bool usePinned = false)
{
	if(!usePinned)
	{
		::free(var); 
	}
	else
	{
		cudaThreadSynchronize();	 // To prevent getting confused over other errors
		cudaError_t cudaMemStatus = cudaFreeHost(var);
		if(cudaMemStatus!=cudaSuccess) ERROR(cudaGetErrorString(cudaMemStatus));
	}
}

// Typesafe re-allocator (convenience)
template<typename T>
T* hostAlloc(T* var, int nelem, bool usePinned = false)
{
	if(!usePinned)
	{
		T* tmp = (T*)realloc(var, nelem*sizeof(T));
		if(tmp == NULL) ERROR("Out of host memory.");
		return tmp;
	}
	else
	{
		cudaThreadSynchronize();   // To prevent getting confused over other errors
		if(var!=NULL) hostFree(var);
		cudaError_t cudaMemStatus = cudaMallocHost((void**)&var,nelem*sizeof(T));
		if(cudaMemStatus!=cudaSuccess) ERROR(cudaGetErrorString(cudaMemStatus));
		return var;
	}
}

//
// Utilities
//

/// trim whitespaces from the beginning and the end of a string
void trim(std::string& str);

// NOTE: The ifdef here is a workaround for CUDA 2.2 device emulation mode bug, where C++
// is disabled in emulation mode. If you want to use the functions below, use them only
// from units compiled with g++
#ifndef __CUDACC__

	#include <valarray>
	void calc_total_energy(const cpu_ensemble &ens, std::valarray<double> &E);

#endif

#endif
