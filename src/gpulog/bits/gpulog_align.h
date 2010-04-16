#ifndef bits_gpulog_align_h__
#define bits_gpulog_align_h__

//
// Struct alignment is handled differently between the CUDA compiler and other
// compilers (e.g. GCC, MS Visual C++ .NET)
//
#ifndef ALIGN
	#define LOCALLY_DEFINED_ALIGN

	#ifdef __CUDACC__
		#define ALIGN(x)  __align__(x)
	#else
		#if defined(_MSC_VER) && (_MSC_VER >= 1300)
			// Visual C++ .NET and later
			#define ALIGN(x) __declspec(align(x))
		#else
			#if defined(__GNUC__)
				// GCC
				#define ALIGN(x)  __attribute__ ((aligned (x)))
			#else
			// all other compilers  
				#define ALIGN(x)
			#endif
		#endif
	#endif
#endif


#endif // bits_gpulog_align_h__
