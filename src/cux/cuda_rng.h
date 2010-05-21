/***************************************************************************
 *   Copyright (C) 2009-2010 by Mario Juric                                *
 *   mjuric@astro.Princeton.EDU                                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef cuda_rng_h__
#define cuda_rng_h__

/****************************************************************************/
/*                                                                          */
/*  A set of multithreaded GPU integrators.                                 */
/*                                                                          */
/****************************************************************************/

/*

Host interface:
	typedef prngs::gpu::mwc cuda_rng;
	cuda_rng rng = cuda_rng::create(seed, nstreams); // nstreams must be >= nthreads
	kernel<<...>>(rng);
	rng.free()

Device interface:
	rng.load();
	...
	u = rng.uniform();
	...
	rng.store();
*/

#include <assert.h>
#include <cmath>

namespace prngs
{
	template<uint32_t state_dim, bool on_gpu>
	struct rng_base
	{
		static const int statewidth = state_dim;
		static int state_width() { return state_dim; }
		static int state_bytes() { return state_dim * sizeof(uint32_t); }

		uint32_t *gstate;
		uint32_t nstreams;

		__host__ void construct_base() { gstate = NULL; nstreams = 0; }

//		#ifdef __CUDACC__
		__device__ void load(uint32_t tid) const
		{
			// load rng state to shared memory
			for(uint32_t i=0; i != statewidth; i++)
			{
				shmem(uint32_t)[i*blockDim.x + threadIdx.x] = gstate[i*nstreams + tid];
/*#if __DEVICE_EMULATION__
				if(threadIdx.x == 0) {
					printf("Thread 0, state[%d]=%u\n", i, shmem(uint32_t)[i*blockDim.x + threadIdx.x]);
				}
#endif*/
			}
		}

		__device__ void store(uint32_t tid) const
		{
			// store rng state to global memory
			for(uint32_t i=0; i != statewidth; i++)
			{
/*			#if __DEVICE_EMULATION__
				printf("tid=%u g=%u s=%u\n", tid, gstate[i*nstreams + tid], shmem(uint32_t)[i*blockDim.x + threadIdx.x]);
			#endif*/
				gstate[i*nstreams + tid] = shmem(uint32_t)[i*blockDim.x + threadIdx.x];
			}
		}
//		#endif

	public:
		// upload the state vector to the GPU
		void upload(const uint32_t *states, const uint32_t nstreams)
		{
			// copy rng stream state to GPU
			free();
			this->nstreams = nstreams;
			if(on_gpu)
			{
				cuxErrCheck( cudaMalloc((void**)&gstate, sizeof(uint32_t)*nstreams*statewidth) );
				cuxErrCheck( cudaMemcpy(gstate, states, sizeof(uint32_t)*nstreams*statewidth, cudaMemcpyHostToDevice) );
			} else {
				gstate = new uint32_t[nstreams*statewidth];
				memcpy(gstate, states, sizeof(uint32_t)*nstreams*statewidth);
			}
		}

		// download the states. The caller must ensure there's enough room in
		// the output buffer.
		int download(uint32_t *states)
		{
			if(!gstate) { return 0; }
			if(on_gpu)
			{
				cuxErrCheck( cudaMemcpy(states, gstate, sizeof(uint32_t)*nstreams*statewidth, cudaMemcpyDeviceToHost) );
			}
			else
			{
				memcpy(states, gstate, sizeof(uint32_t)*nstreams*statewidth);
			}
			return nstreams*statewidth;
		}

		// Rounds up v to nearest integer divisable by mod
		uint32_t roundUpModulo(uint32_t v, uint32_t mod)
		{
			uint32_t r = v % mod;
			uint32_t pitch = r ? v + (mod-r) : v;
			return pitch;
		}

	public:
		// free the GPU state vector
		void free()
		{
			if(!gstate) { return; }

			if(on_gpu)
			{
				cuxErrCheck( cudaFree(gstate) );
			}
			else
			{
				delete [] gstate;
			}
			gstate = NULL;
		}

		void srand(uint32_t seed, uint32_t nstreams_=0)
		{
			if(nstreams_) { nstreams = roundUpModulo(nstreams_, 64); }
			uint32_t *streams = new uint32_t[nstreams*statewidth];

			srand48(seed);
			for(uint32_t i = 0; i != nstreams; i++)
			{
				streams[i] = lrand48();
			}
			upload(streams, nstreams);
			delete [] streams;
		}
		static const char *name() { return "unknown"; }

	};

	template<typename rng_impl>
	struct rng : public rng_impl
	{
	public:
		static rng create(uint32_t seed, uint32_t nstreams) { rng r; r.construct_base(); r.srand(seed, nstreams); return r; }
		static rng create() { rng r; r.construct_base(); return r; }

		//typedef rng_noctor<rng_impl> constant;
		typedef rng constant;

//		#ifdef __CUDACC__
		__device__ float uniform_pos() const
		{
			float x;
			do { x = this->uniform(); } while (x == 0.f);
			return x;
		}

		__device__ float gaussian(const float sigma) const
		{
			float x, y, r2;

			do
			{
				/* choose x,y in uniform square (-1,-1) to (+1,+1) */
				x = -1.f + 2.f * this->uniform();
				y = -1.f + 2.f * this->uniform();

				/* see if it is in the unit circle */
				r2 = x * x + y * y;
			}
			while (r2 > 1.0f || r2 == 0.f || r2 == 4.f);

			/* Box-Muller transform */
			return sigma * y * sqrt (-2.0f * logf (r2) / r2);
		}

		// adapted from GNU Scientific Library (GSL)
		__device__ float gamma_large(const float a) const
		{
			/* Works only if a > 1, and is most efficient if a is large
		
				 This algorithm, reported in Knuth, is attributed to Ahrens.	A
				 faster one, we are told, can be found in: J. H. Ahrens and
				 U. Dieter, Computing 12 (1974) 223-246.	*/

			float sqa, x, y, v;
			sqa = sqrtf (2 * a - 1);
			do
			{
				do
				{
					y = tanf (M_PI * this->uniform());
					x = sqa * y + a - 1;
				}
				while (x <= 0);
				v = this->uniform();
			}
			while(v > (1 + y*y) * expf ((a - 1) * logf (x / (a - 1)) - sqa * y));

			return x;
		}

		// adapted from GNU Scientific Library (GSL)
		__device__ float gamma_int(const unsigned int a) const
		{
			if (a < 3)
			{
				unsigned int i;
				float prod = 1;

				for (i = 0; i < a; i++)
				{
					prod *= this->uniform_pos();
				}
		
				/* Note: for 12 iterations we are safe against underflow, since
				 the smallest positive random number is O(2^-32). This means
				 the smallest possible product is 2^(-12*32) = 10^-116 which
				 is within the range of double precision. */
				/* mjuric: I modified the number of iterations (3) to reflect
				 the change from double to single precision arithmetic. */
				return -logf(prod);
			}
			else
			{
				return gamma_large((float)a);
			}
		}

		/* Taken from GSL */
		__device__ float gamma_frac(const float a) const
		{
			/* This is exercise 16 from Knuth; see page 135, and the solution is
				 on page 551.	*/
		
			float p, q, x, u, v;
			p = M_E / (a + M_E);
			do
			{
				u = this->uniform();
				v = uniform_pos();

				if (u < p)
				{
					x = exp ((1 / a) * log (v));
					q = exp (-x);
				}
				else
				{
					x = 1 - log (v);
					q = exp ((a - 1) * log (x));
				}
			}
			while (this->uniform() >= q);

			return x;
		}

		/* The Gamma distribution of order a>0 is defined by:

		   p(x) dx = {1 / \Gamma(a) b^a } x^{a-1} e^{-x/b} dx

		   for x>0.  If X and Y are independent gamma-distributed random
		   variables of order a1 and a2 with the same scale parameter b, then
		   X+Y has gamma distribution of order a1+a2.

		   The algorithms below are from Knuth, vol 2, 2nd ed, p. 129. */

		__device__ float gamma(const float a, const float b) const
		{
			/* assume a > 0 */
			unsigned int na = floor (a);

			if(a >= UINT_MAX)
			{
				return b * (gamma_large (floor (a)) + gamma_frac(a - floor(a)));
			}
			else if (a == na)
			{
				return b * gamma_int(na);
			}
			else if (na == 0)
			{
				return b * gamma_frac(a);
			}
			else
			{
				return b * (gamma_int(na) + gamma_frac(a - na)) ;
			}
		}

		__device__ float beta(const float a, const float b) const
		{
			float x1 = gamma(a, 1.0);
			float x2 = gamma(b, 1.0);

			return x1 / (x1 + x2);
		}

		/* Adapted from binomial_knuth from GSL */
		__device__ unsigned int binomial(float p, unsigned int n) const
		{
			unsigned int i, a, b, k = 0;

			while (n > 10)				/* This parameter is tunable */
			{
				float X;
				a = 1 + (n / 2);
				b = 1 + n - a;
		
				X = beta ((float) a, (float) b);
		
				if (X >= p)
				{
					n = a - 1;
					p /= X;
				}
				else
				{
					k += a;
					n = b - 1;
					p = (p - X) / (1 - X);
				}
			}
		
			for (i = 0; i < n; i++)
			{
				float u = this->uniform();
				if (u < p)
					k++;
			}
			return k;
		}


		// Adapted from GNU Scientific Library (GSL)
		__device__ unsigned int poisson(float mu) const
		{
			float emu;
			float prod = 1.0;
			unsigned int k = 0;
#if 1
			while (mu > 10)
			{
				unsigned int m = mu * (7.0 / 8.0);
				float X = gamma_int(m);
		
				if (X >= mu)
				{
					return k + binomial(mu / X, m - 1);
				}
				else
				{
					k += m;
					mu -= X; 
				}
			}
#endif
			/* This following method works well when mu is small */
			emu = expf(-mu);
			do
			{
				prod *= this->uniform();
				k++;
			}
			while(prod>emu);

			return k - 1;
		}
//		#endif
	};

	template<bool on_gpu>
	struct ran0_impl : public rng_base<1, on_gpu>
	{
		static const int IA = 16807;
		static const int IM = 2147483647;
		static const float AM = (1.0f/IM);
		static const int IQ = 127773;
		static const int IR = 2836;
		static const int MASK = 123459876;

		// An ultra-simple random number generator (straight out of NR)
//		#ifdef __CUDACC__
		__device__ float uniform() const
		{
			int idum = shmem(uint32_t)[threadIdx.x];

			int k;
			float ans;

			idum ^= MASK;			// XORing with MASK allows use of zero and other
			k=idum/IQ;			// simple bit patterns for idum.
			idum=IA*(idum-k*IQ)-IR*k;	// Compute idum=(IA*idum) % IM without over-
			if (idum < 0) idum += IM;	// flows by Schrage's method.
			ans=AM*idum; 			// Convert idum to a floating result.
			idum ^= MASK; 			// Unmask before return.

			shmem(uint32_t)[threadIdx.x] = idum;
			return ans;
		}
//		#endif

		static const char *name() { return "ran0"; }
	};

	template<bool on_gpu>
	struct mwc_impl : public rng_base<3, on_gpu>
	{
//		#ifdef __CUDACC__
#if 0	// if this is on, the integrator state is stored in registers
		uint32_t a, c, xn;
		__device__ void load(uint32_t tid)
		{
			a  = gstate[tid];
			c  = gstate[nstreams + tid];
			xn = gstate[2*nstreams + tid];
		}

		__device__ void store(uint32_t tid)
		{
			gstate[tid] = a;
			gstate[nstreams + tid] = c;
			gstate[2*nstreams + tid] = xn;
		}
#else
		#define a  (shmem(uint32_t)[               threadIdx.x])
		#define c  (shmem(uint32_t)[  blockDim.x + threadIdx.x])
		#define xn (shmem(uint32_t)[2*blockDim.x + threadIdx.x])
#endif
		// An ultra-simple random number generator (straight out of NR)
		__device__ float uniform() const
		{

			/*
				Marsaglia's Multiply-With-Carry RNG. For theory and details see:
			
					http://www.stat.fsu.edu/pub/diehard/cdrom/pscript/mwc1.ps
					http://www.ms.uky.edu/~mai/RandomNumber
					http://www.ast.cam.ac.uk/~stg20/cuda/random/index.html
			*/
/*			#if __DEVICE_EMULATION__
				printf("ti.x=%d %u %u %u\n", threadIdx.x, a, c, xn);
			#endif*/
			uint64_t xnew = (uint64_t)a*xn + c;
			c = xnew >> 32;
			xn = (xnew << 32) >> 32;
/*			#if __DEVICE_EMULATION__
				printf("ti.x=%d %u %u %u\n", threadIdx.x, a, c, xn);
			#endif*/
			return 2.32830643708e-10f * xn;

			#undef a
			#undef c
			#undef xn
		}
//		#endif

		void srand(uint32_t seed, uint32_t nstreams=0, const char *safeprimes_file = NULL)
		{
			if(nstreams) { nstreams = roundUpModulo(nstreams, 64); }
			uint32_t *streams = new uint32_t[3*nstreams];

			assert(nstreams <= 1<<16);

			// initialize CPU streams
			if(safeprimes_file == NULL) { safeprimes_file = "safeprimes32.txt"; }
			FILE *f = fopen(safeprimes_file, "r");
			assert(f != NULL);
			srand48(seed);
			//printf("seed=%u\n", seed);
			for(uint32_t i = 0; i != nstreams; i++)
			{
				unsigned int prime;
				unsigned long d1, d2;

				fscanf(f, "%u %lu %lu\n", &prime, &d1, &d2);
				streams[i] = prime;					// multiplier
				streams[  nstreams + i] = (int)(drand48() * prime);	// initial carry (nas to be < multiplier)
				streams[2*nstreams + i] = lrand48();			// initial x
			}
			fclose(f);

			this->upload(streams, nstreams);
			delete [] streams;
		}

		static const char *name() { return "mwc"; }
	};

	template<bool on_gpu>
	struct taus2_impl : public rng_base<3, on_gpu>
	{
		#ifdef __CUDACC__
		__device__ float uniform() const
		{
			/*
				taus2 generator, from:
				http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html
			*/
			#define s1     (shmem(uint32_t)[               threadIdx.x])
			#define s2     (shmem(uint32_t)[  blockDim.x + threadIdx.x])
			#define s3     (shmem(uint32_t)[2*blockDim.x + threadIdx.x])

			s1 = (((s1 & 4294967294)<<12)^(((s1<<13)^s1)>>19));
			s2 = (((s2 & 4294967288)<< 4)^(((s2<< 2)^s2)>>25));
			s3 = (((s3 & 4294967280)<<17)^(((s3<< 3)^s3)>>11));
			uint32_t xn = (s1 ^ s2 ^ s3);

			return 2.32830643708e-10f * xn;

			#undef s1
			#undef s2
			#undef s3
		}
		#endif

		static const char *name() { return "taus2"; }
	};

	template<bool on_gpu>
	struct rand48_impl : public rng_base<3, on_gpu>
	{
		// Adapted from:
		// 	http://forums.nvidia.com/index.php?act=attach&type=post&id=9512
		//

		#define statex  (shmem(uint32_t)[               threadIdx.x])
		#define statey  (shmem(uint32_t)[  blockDim.x + threadIdx.x])

		static const uint32_t Ax = 0x00ECE66D; // 0x5DE ECE66D = 25214903917;
		static const uint32_t Ay = 0x000005DE; // 0x5DE ECE66D = 25214903917;
		static const uint32_t Cx = 0x0000000B; // 0xB = 11
		static const uint32_t Cy = 0x00000000; // 0xB = 11

		#ifdef __CUDACC__
		__device__ inline void rand48_iterate() const
		{

			// state0 is 2x 24bit to handle overflows optimally, i.e.
			// in one operation.

			// the multiplication commands however give the low and hi 32 bit,
			// which have to be converted as follows:
			// 48bit in bytes = ABCD EF (space marks 32bit boundary)
			// R0             = ABC
			// R1             =    D EF

			unsigned int R0, R1;
	
			// low 24-bit multiplication
			const unsigned int lo00 = __umul24(statex, Ax);
			const unsigned int hi00 = __umulhi(statex, Ax);

			// 24bit distribution of 32bit multiplication results
			R0 = (lo00 & 0xFFFFFF);
			R1 = (lo00 >> 24) | (hi00 << 8);
	
			R0 += Cx; R1 += Cy;
	
			// transfer overflows
			R1 += (R0 >> 24);
			R0 &= 0xFFFFFF;
	
			// cross-terms, low/hi 24-bit multiplication
			R1 += __umul24(statey, Ax);
			R1 += __umul24(statex, Ay);

			R1 &= 0xFFFFFF;

			statex = R0;
			statey = R1;
		}

		__device__ inline int rand48_nextInt() const
		{
			// get upper 31 (!) bits of the 2x 24bits
			int res = ( statex >> 17 ) | ( statey << 7 );
			rand48_iterate();
			return res;
		}

		// returns a float in the range [0, 1)
		__device__ inline float uniform() const
		{
			// use only upper 24 bits since floating point has 24 bit significand
			// (ref: Java random documentation)
			float res = (float)statey / (float)(1<<24);
			rand48_iterate();
			return res;
		}
		#endif

		#undef statex
		#undef statey

		void srand(uint32_t seed, int nstreams=0)
		{
			if(nstreams) { nstreams = roundUpModulo(nstreams, 64); }
			uint32_t *streams = new uint32_t[this->statewidth*nstreams];

			for(int i = 0; i != nstreams; i++)
			{
				unsigned long long x = (((unsigned long long)(seed+i)) << 16) | 0x330E;
				streams[i]           =         x & 0xFFFFFFLL;
				streams[nstreams+i]  = (x >> 24) & 0xFFFFFFLL;
			}

			this->upload(streams, nstreams);
			delete [] streams;
		}

		static const char *name() { return "rand48"; }
	};
	
	namespace gpu
	{
		typedef rng<ran0_impl<true> >   ran0;
		typedef rng<mwc_impl<true> >    mwc;
		typedef rng<taus2_impl<true> >  taus2;
		typedef rng<rand48_impl<true> > rand48;
	}
	
	namespace cpu
	{
		typedef rng<ran0_impl<false> >   ran0;
		typedef rng<mwc_impl<false> >    mwc;
		typedef rng<taus2_impl<false> >  taus2;
		typedef rng<rand48_impl<false> > rand48;
	}
}

#endif // cuda_rng_h__
