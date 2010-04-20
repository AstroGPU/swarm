#ifndef bits_gpulog_macro_cleanup_h__
#define bits_gpulog_macro_cleanup_h__

	// Macro cleanup
	#ifdef LOCALLY_DEFINED_ALIGN
		#undef ALIGN
		#undef LOCALLY_DEFINED_ALIGN
	#endif

	#undef ISARRAY
//	#undef SCALAR
	#undef ISUNSPEC
	#undef A
	#undef M
	#undef XSTART
	#undef XSIZEOF
	#undef ASTART
	#undef SIZEOF
	#undef ADDR
//	#undef PTR_T

#endif // bits_gpulog_macro_cleanup_h__
