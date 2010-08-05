//
// !!! MACHINE GENERATED FILE; DO NOT EDIT !!!
// Generated using ./gen_gpulog_write.pl printf
//
	template<typename L, int N>
#ifdef __CUDACC__
        __device__
#endif
		inline void lprintf(L &log, const char (&fmt)[N])
	{
		log.write(gpulog::MSG_PRINTF, fmt);
	}
	template<typename L, int N, typename T1>
#ifdef __CUDACC__
        __device__
#endif
		inline void lprintf(L &log, const char (&fmt)[N], const T1 &v1)
	{
		log.write(gpulog::MSG_PRINTF, fmt, f2d(T1, v1));
	}
	template<typename L, int N, typename T1, typename T2>
#ifdef __CUDACC__
        __device__
#endif
		inline void lprintf(L &log, const char (&fmt)[N], const T1 &v1, const T2 &v2)
	{
		log.write(gpulog::MSG_PRINTF, fmt, f2d(T1, v1), f2d(T2, v2));
	}
	template<typename L, int N, typename T1, typename T2, typename T3>
#ifdef __CUDACC__
        __device__
#endif
		inline void lprintf(L &log, const char (&fmt)[N], const T1 &v1, const T2 &v2, const T3 &v3)
	{
		log.write(gpulog::MSG_PRINTF, fmt, f2d(T1, v1), f2d(T2, v2), f2d(T3, v3));
	}
	template<typename L, int N, typename T1, typename T2, typename T3, typename T4>
#ifdef __CUDACC__
        __device__
#endif
		inline void lprintf(L &log, const char (&fmt)[N], const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4)
	{
		log.write(gpulog::MSG_PRINTF, fmt, f2d(T1, v1), f2d(T2, v2), f2d(T3, v3), f2d(T4, v4));
	}
	template<typename L, int N, typename T1, typename T2, typename T3, typename T4, typename T5>
#ifdef __CUDACC__
        __device__
#endif
		inline void lprintf(L &log, const char (&fmt)[N], const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5)
	{
		log.write(gpulog::MSG_PRINTF, fmt, f2d(T1, v1), f2d(T2, v2), f2d(T3, v3), f2d(T4, v4), f2d(T5, v5));
	}
	template<typename L, int N, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
#ifdef __CUDACC__
        __device__
#endif
		inline void lprintf(L &log, const char (&fmt)[N], const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6)
	{
		log.write(gpulog::MSG_PRINTF, fmt, f2d(T1, v1), f2d(T2, v2), f2d(T3, v3), f2d(T4, v4), f2d(T5, v5), f2d(T6, v6));
	}
	template<typename L, int N, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
#ifdef __CUDACC__
        __device__
#endif
		inline void lprintf(L &log, const char (&fmt)[N], const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7)
	{
		log.write(gpulog::MSG_PRINTF, fmt, f2d(T1, v1), f2d(T2, v2), f2d(T3, v3), f2d(T4, v4), f2d(T5, v5), f2d(T6, v6), f2d(T7, v7));
	}
	template<typename L, int N, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
#ifdef __CUDACC__
        __device__
#endif
		inline void lprintf(L &log, const char (&fmt)[N], const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8)
	{
		log.write(gpulog::MSG_PRINTF, fmt, f2d(T1, v1), f2d(T2, v2), f2d(T3, v3), f2d(T4, v4), f2d(T5, v5), f2d(T6, v6), f2d(T7, v7), f2d(T8, v8));
	}
	template<typename L, int N, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
#ifdef __CUDACC__
        __device__
#endif
		inline void lprintf(L &log, const char (&fmt)[N], const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8, const T9 &v9)
	{
		log.write(gpulog::MSG_PRINTF, fmt, f2d(T1, v1), f2d(T2, v2), f2d(T3, v3), f2d(T4, v4), f2d(T5, v5), f2d(T6, v6), f2d(T7, v7), f2d(T8, v8), f2d(T9, v9));
	}
