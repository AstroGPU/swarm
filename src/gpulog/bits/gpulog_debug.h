#ifndef bits_gpulog_debug_h__
#define bits_gpulog_debug_h__

//
// Debugging macros. Define GPULOG_DEBUG=1 to turn on
//
#if GPULOG_DEBUG
	#define DBG(x) x
#else
	#define DBG(x)
#endif

#if !__CUDACC__
	#include <iostream>
	#include <cassert>
	#include <sstream>

	#define DHOST(x) DBG(x)
	#define DGPU(x)
	
	#define dev_assert(x) assert(x)
#else
	#define DHOST(x)	/* This hides the stuff that nvcc can't compile */

	#if __DEVICE_EMULATION__
		#define DGPU(x) DBG(x)
		// #define dev_assert(x) assert(x) /* has to be disabled on CUDA 2.3 or compilation fails */
		#define dev_assert(x) { if(!(x)) { printf("Assertion failed: " #x "\n"); exit(-1); } } /* This is a CUDA 2.3 compatible replacement for assert(); */
		// #define dev_assert(x)
	#else
		#define DGPU(x)
		#define dev_assert(x)
	#endif
#endif

#endif // bits_gpulog_debug_h
