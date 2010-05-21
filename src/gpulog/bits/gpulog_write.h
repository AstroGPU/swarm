//
// !!! MACHINE GENERATED FILE; DO NOT EDIT !!!
// Generated using ../scripts/gen_gpulog_write.pl write
//
template<typename T1>
	__host__ __device__ inline PTR_T(SCALAR(T1)) write(const int recid, const T1 &v1)
	{
		typedef internal::pktsize<header, T1> P;
		//P::dump();

		// allocate and test for end-of-buffer
		int len = P::len_with_padding(v1);
		int at = A::atomicAdd(this->at, len);
		if(has_overflowed(at + len)) 
		{
			A::atomicAdd(this->at, -len);
			return NULL;
		}
		char *ptr = buffer + at;

		// write
		header v0(recid, len);
		P::IO0::put(ptr, v0, P::begin0, P::len0);
		P::IO1::put(ptr, v1, P::begin1, P::len1);

#if ARGINFO
		P::store_arginfo(ptr, v1);
#endif

		DHOST( std::cerr << "Total packet len = " << len << "\n"; )
		return (SCALAR(T1)*)(ptr + P::begin1);
	}

template<typename T1, typename T2>
	__host__ __device__ inline PTR_T(SCALAR(T2)) write(const int recid, const T1 &v1, const T2 &v2)
	{
		typedef internal::pktsize<header, T1, T2> P;
		//P::dump();

		// allocate and test for end-of-buffer
		int len = P::len_with_padding(v2);
		int at = A::atomicAdd(this->at, len);
		if(has_overflowed(at + len)) 
		{
			A::atomicAdd(this->at, -len);
			return NULL;
		}
		char *ptr = buffer + at;

		// write
		header v0(recid, len);
		P::IO0::put(ptr, v0, P::begin0, P::len0);
		P::IO1::put(ptr, v1, P::begin1, P::len1);
		P::IO2::put(ptr, v2, P::begin2, P::len2);

#if ARGINFO
		P::store_arginfo(ptr, v2);
#endif

		DHOST( std::cerr << "Total packet len = " << len << "\n"; )
		return (SCALAR(T2)*)(ptr + P::begin2);
	}

template<typename T1, typename T2, typename T3>
	__host__ __device__ inline PTR_T(SCALAR(T3)) write(const int recid, const T1 &v1, const T2 &v2, const T3 &v3)
	{
		typedef internal::pktsize<header, T1, T2, T3> P;
		//P::dump();

		// allocate and test for end-of-buffer
		int len = P::len_with_padding(v3);
		int at = A::atomicAdd(this->at, len);
		if(has_overflowed(at + len)) 
		{
			A::atomicAdd(this->at, -len);
			return NULL;
		}
		char *ptr = buffer + at;

		// write
		header v0(recid, len);
		P::IO0::put(ptr, v0, P::begin0, P::len0);
		P::IO1::put(ptr, v1, P::begin1, P::len1);
		P::IO2::put(ptr, v2, P::begin2, P::len2);
		P::IO3::put(ptr, v3, P::begin3, P::len3);

#if ARGINFO
		P::store_arginfo(ptr, v3);
#endif

		DHOST( std::cerr << "Total packet len = " << len << "\n"; )
		return (SCALAR(T3)*)(ptr + P::begin3);
	}

template<typename T1, typename T2, typename T3, typename T4>
	__host__ __device__ inline PTR_T(SCALAR(T4)) write(const int recid, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4)
	{
		typedef internal::pktsize<header, T1, T2, T3, T4> P;
		//P::dump();

		// allocate and test for end-of-buffer
		int len = P::len_with_padding(v4);
		int at = A::atomicAdd(this->at, len);
		if(has_overflowed(at + len)) 
		{
			A::atomicAdd(this->at, -len);
			return NULL;
		}
		char *ptr = buffer + at;

		// write
		header v0(recid, len);
		P::IO0::put(ptr, v0, P::begin0, P::len0);
		P::IO1::put(ptr, v1, P::begin1, P::len1);
		P::IO2::put(ptr, v2, P::begin2, P::len2);
		P::IO3::put(ptr, v3, P::begin3, P::len3);
		P::IO4::put(ptr, v4, P::begin4, P::len4);

#if ARGINFO
		P::store_arginfo(ptr, v4);
#endif

		DHOST( std::cerr << "Total packet len = " << len << "\n"; )
		return (SCALAR(T4)*)(ptr + P::begin4);
	}

template<typename T1, typename T2, typename T3, typename T4, typename T5>
	__host__ __device__ inline PTR_T(SCALAR(T5)) write(const int recid, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5)
	{
		typedef internal::pktsize<header, T1, T2, T3, T4, T5> P;
		//P::dump();

		// allocate and test for end-of-buffer
		int len = P::len_with_padding(v5);
		int at = A::atomicAdd(this->at, len);
		if(has_overflowed(at + len)) 
		{
			A::atomicAdd(this->at, -len);
			return NULL;
		}
		char *ptr = buffer + at;

		// write
		header v0(recid, len);
		P::IO0::put(ptr, v0, P::begin0, P::len0);
		P::IO1::put(ptr, v1, P::begin1, P::len1);
		P::IO2::put(ptr, v2, P::begin2, P::len2);
		P::IO3::put(ptr, v3, P::begin3, P::len3);
		P::IO4::put(ptr, v4, P::begin4, P::len4);
		P::IO5::put(ptr, v5, P::begin5, P::len5);

#if ARGINFO
		P::store_arginfo(ptr, v5);
#endif

		DHOST( std::cerr << "Total packet len = " << len << "\n"; )
		return (SCALAR(T5)*)(ptr + P::begin5);
	}

template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
	__host__ __device__ inline PTR_T(SCALAR(T6)) write(const int recid, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6)
	{
		typedef internal::pktsize<header, T1, T2, T3, T4, T5, T6> P;
		//P::dump();

		// allocate and test for end-of-buffer
		int len = P::len_with_padding(v6);
		int at = A::atomicAdd(this->at, len);
		if(has_overflowed(at + len)) 
		{
			A::atomicAdd(this->at, -len);
			return NULL;
		}
		char *ptr = buffer + at;

		// write
		header v0(recid, len);
		P::IO0::put(ptr, v0, P::begin0, P::len0);
		P::IO1::put(ptr, v1, P::begin1, P::len1);
		P::IO2::put(ptr, v2, P::begin2, P::len2);
		P::IO3::put(ptr, v3, P::begin3, P::len3);
		P::IO4::put(ptr, v4, P::begin4, P::len4);
		P::IO5::put(ptr, v5, P::begin5, P::len5);
		P::IO6::put(ptr, v6, P::begin6, P::len6);

#if ARGINFO
		P::store_arginfo(ptr, v6);
#endif

		DHOST( std::cerr << "Total packet len = " << len << "\n"; )
		return (SCALAR(T6)*)(ptr + P::begin6);
	}

template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
	__host__ __device__ inline PTR_T(SCALAR(T7)) write(const int recid, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7)
	{
		typedef internal::pktsize<header, T1, T2, T3, T4, T5, T6, T7> P;
		//P::dump();

		// allocate and test for end-of-buffer
		int len = P::len_with_padding(v7);
		int at = A::atomicAdd(this->at, len);
		if(has_overflowed(at + len)) 
		{
			A::atomicAdd(this->at, -len);
			return NULL;
		}
		char *ptr = buffer + at;

		// write
		header v0(recid, len);
		P::IO0::put(ptr, v0, P::begin0, P::len0);
		P::IO1::put(ptr, v1, P::begin1, P::len1);
		P::IO2::put(ptr, v2, P::begin2, P::len2);
		P::IO3::put(ptr, v3, P::begin3, P::len3);
		P::IO4::put(ptr, v4, P::begin4, P::len4);
		P::IO5::put(ptr, v5, P::begin5, P::len5);
		P::IO6::put(ptr, v6, P::begin6, P::len6);
		P::IO7::put(ptr, v7, P::begin7, P::len7);

#if ARGINFO
		P::store_arginfo(ptr, v7);
#endif

		DHOST( std::cerr << "Total packet len = " << len << "\n"; )
		return (SCALAR(T7)*)(ptr + P::begin7);
	}

template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
	__host__ __device__ inline PTR_T(SCALAR(T8)) write(const int recid, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8)
	{
		typedef internal::pktsize<header, T1, T2, T3, T4, T5, T6, T7, T8> P;
		//P::dump();

		// allocate and test for end-of-buffer
		int len = P::len_with_padding(v8);
		int at = A::atomicAdd(this->at, len);
		if(has_overflowed(at + len)) 
		{
			A::atomicAdd(this->at, -len);
			return NULL;
		}
		char *ptr = buffer + at;

		// write
		header v0(recid, len);
		P::IO0::put(ptr, v0, P::begin0, P::len0);
		P::IO1::put(ptr, v1, P::begin1, P::len1);
		P::IO2::put(ptr, v2, P::begin2, P::len2);
		P::IO3::put(ptr, v3, P::begin3, P::len3);
		P::IO4::put(ptr, v4, P::begin4, P::len4);
		P::IO5::put(ptr, v5, P::begin5, P::len5);
		P::IO6::put(ptr, v6, P::begin6, P::len6);
		P::IO7::put(ptr, v7, P::begin7, P::len7);
		P::IO8::put(ptr, v8, P::begin8, P::len8);

#if ARGINFO
		P::store_arginfo(ptr, v8);
#endif

		DHOST( std::cerr << "Total packet len = " << len << "\n"; )
		return (SCALAR(T8)*)(ptr + P::begin8);
	}

template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
	__host__ __device__ inline PTR_T(SCALAR(T9)) write(const int recid, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8, const T9 &v9)
	{
		typedef internal::pktsize<header, T1, T2, T3, T4, T5, T6, T7, T8, T9> P;
		//P::dump();

		// allocate and test for end-of-buffer
		int len = P::len_with_padding(v9);
		int at = A::atomicAdd(this->at, len);
		if(has_overflowed(at + len)) 
		{
			A::atomicAdd(this->at, -len);
			return NULL;
		}
		char *ptr = buffer + at;

		// write
		header v0(recid, len);
		P::IO0::put(ptr, v0, P::begin0, P::len0);
		P::IO1::put(ptr, v1, P::begin1, P::len1);
		P::IO2::put(ptr, v2, P::begin2, P::len2);
		P::IO3::put(ptr, v3, P::begin3, P::len3);
		P::IO4::put(ptr, v4, P::begin4, P::len4);
		P::IO5::put(ptr, v5, P::begin5, P::len5);
		P::IO6::put(ptr, v6, P::begin6, P::len6);
		P::IO7::put(ptr, v7, P::begin7, P::len7);
		P::IO8::put(ptr, v8, P::begin8, P::len8);
		P::IO9::put(ptr, v9, P::begin9, P::len9);

#if ARGINFO
		P::store_arginfo(ptr, v9);
#endif

		DHOST( std::cerr << "Total packet len = " << len << "\n"; )
		return (SCALAR(T9)*)(ptr + P::begin9);
	}

template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
	__host__ __device__ inline PTR_T(SCALAR(T10)) write(const int recid, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8, const T9 &v9, const T10 &v10)
	{
		typedef internal::pktsize<header, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10> P;
		//P::dump();

		// allocate and test for end-of-buffer
		int len = P::len_with_padding(v10);
		int at = A::atomicAdd(this->at, len);
		if(has_overflowed(at + len)) 
		{
			A::atomicAdd(this->at, -len);
			return NULL;
		}
		char *ptr = buffer + at;

		// write
		header v0(recid, len);
		P::IO0::put(ptr, v0, P::begin0, P::len0);
		P::IO1::put(ptr, v1, P::begin1, P::len1);
		P::IO2::put(ptr, v2, P::begin2, P::len2);
		P::IO3::put(ptr, v3, P::begin3, P::len3);
		P::IO4::put(ptr, v4, P::begin4, P::len4);
		P::IO5::put(ptr, v5, P::begin5, P::len5);
		P::IO6::put(ptr, v6, P::begin6, P::len6);
		P::IO7::put(ptr, v7, P::begin7, P::len7);
		P::IO8::put(ptr, v8, P::begin8, P::len8);
		P::IO9::put(ptr, v9, P::begin9, P::len9);
		P::IO10::put(ptr, v10, P::begin10, P::len10);

#if ARGINFO
		P::store_arginfo(ptr, v10);
#endif

		DHOST( std::cerr << "Total packet len = " << len << "\n"; )
		return (SCALAR(T10)*)(ptr + P::begin10);
	}

