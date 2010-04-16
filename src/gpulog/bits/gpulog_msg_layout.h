#ifndef bits_gpulog_msg_layout_h__
#define bits_gpulog_msg_layout_h__

/*
	class template pktsize<> and supporting classes and templates that
	enable the compile-time calculation of log record layout (byte
	offsets to which the data will be stored).
*/

namespace gpulog
{
	namespace internal
	{

	#if !__CUDACC__
	//
	// Aux operators for debugging
	//
	template<typename T>
	inline std::ostream &operator <<(std::ostream &out, const array<T> &v) { return out << "n=" << v.nelem; }
	inline std::ostream &operator <<(std::ostream &out, const header &h)   { return out << "msgid=" << h.msgid << " len=" << h.len IFARGINFO( << " nargs=" << h.nargs << " infos=" << h.infos); }
	inline std::ostream &operator <<(std::ostream &out, const Tunspec &tu) { return out; }
	inline std::ostream &operator <<(std::ostream &out, const arginfo &h)
	{
		out << "arg=" << h.arg << " align=" << h.align << " size=" << h.size << " dim=" << h.dim;
		out << " isarr=" << h.isarray << " nelem=" << h.nelem;;
		out << " begin=" << h.begin << " len=" << h.len;
		return out;
	}
	#endif

	//
	// Write support templates
	//	- scalar PODs
	//	- pointers (not allowed, static assertion)
	//	- presized arrays
	//	- array<>s (allocation)
	//

	// i/o support: scalar PODs
	template<typename T> struct argio
	{
		__host__ __device__ static inline void put(char *ptr, const T &x, int start, int datalen)
		{
			DHOST( std::cerr << "Writing [" << x << "] start=" << start << " len=" << datalen << "\n" );
			DGPU( printf("Writing start=%d len=%d\n", start, datalen); );
			*(T*)(ptr + start) = x;
		}
	};

	// force a compiler error if the user attempts to serialize a pointer.
	template<typename T> struct argio<T*>
	{
		__host__ __device__ static inline void put(char *ptr, const T *x, int start, int datalen)
		{
			STATIC_ASSERTION_FAILED__Pointer_serialization_is_not_allowed__Use_fixed_size_C_arrays_or_array_template_instead_____(x);
		}
	};

	// presized array read/write specialization
	template<typename T, int N> struct argio<T[N]>
	{
		// array write specialization
		__host__ __device__ static inline void put(char *ptr, const T *x, int start, int datalen)
		{
			DHOST( std::cerr << "Writing presized array [" << x << "] start=" << start << " len=" << datalen << "\n" );
			DGPU( printf("Writing presized array start=%d len=%d\n", start, datalen); );
			ptr += start;
			for(int i = 0; i != datalen / sizeof(T); i++)
			{
				((T*)ptr)[i] = x[i];
			}
		}
	};

	// i/o support: unbound array (array<>) specialization
	template<typename T> struct argio<array<T> >
	{
		__host__ __device__ static inline void put(char *ptr, const array<T> &x, int start, int datalen)
		{
			// Do nothing. The user must write the array data.
			DHOST( std::cerr << "Allocating array [" << x << "] start=" << start << " element_size=" << datalen << " v[0]= " << *(T*)(ptr + start) << "\n" );
			DGPU( printf("Allocating array start=%d element_size=%d\n", start, datalen); );
		}
	};


	//
	// Compile-time record layout computation machinery
	//

	#define ASTART(at, a)		(at & (a-1) ? (at & ~(a-1)) + a : at)		/* Aligned start address closest but >= than at, for type T */

	// CUDA 2.2/2.3 compilation speedup hack -- otherwise (if ASTART is called directly), nvcc
	// treats static const int vars as _macros_ and attempts to expand them recursively 10+ times (!)
	// The compilation succeeds, but lasts forever (~half a minute)
	template<int at, int a> struct ata { static const int value = ASTART(at, a); };
	#define ASTARTx(at, a) (ata<at, a>::value)

	#define ADDR(beg, k)					\
		static const int align##k = ALIGNOF(T##k);	\
		static const int begin##k = ASTARTx(beg, align##k); \
		static const int   len##k = SIZEOF(T##k);	\
		static const int   end##k = begin##k + len##k;	\
		typedef argio<T##k> IO##k;			\
		__host__ __device__ static inline void dump##k() { if(ISUNSPEC(T##k)) { return; }; DHOST( std::cerr << "arg " << k << ": begin=" << begin##k << " end=" << end##k << " len=" << len##k << " isarray=" << ISARRAY(T##k) << "\n"; ); } \
		__host__ __device__ static inline void get_arginfo##k(arginfo *ai, int nelem) \
		{ \
			typedef ttrait<T##k> TT; \
			if(!ISUNSPEC(T##k)) \
			{  \
				ai->arg = k; \
				ai->align = TT::align; \
				ai->size = TT::size; \
				ai->dim = TT::dim; \
				ai->isarray = TT::isarr; \
				ai->begin = begin##k; \
				ai->nelem = ai->isarray ? nelem : 1; \
				ai->len = len##k * ai->nelem; \
			} \
		}

	/*
		struct template to compile-time compute (properly aligned) offsets and sizes
		of passed types. Used in conjunction with write function templates.
	*/
	template <
		typename T0,
		typename T1 = Tunspec, typename T2 = Tunspec, typename T3 = Tunspec, typename T4 = Tunspec, typename T5 = Tunspec, 
		typename T6 = Tunspec, typename T7 = Tunspec, typename T8 = Tunspec, typename T9 = Tunspec, typename T10 = Tunspec
	>
	struct pktsize
	{
		// number of arguments
		static const int nargs = 
		          !ISUNSPEC(T0) +
			+ !ISUNSPEC(T1) + !ISUNSPEC(T2) + !ISUNSPEC(T3) + !ISUNSPEC(T4) + !ISUNSPEC(T5) + 
			+ !ISUNSPEC(T6) + !ISUNSPEC(T7) + !ISUNSPEC(T8) + !ISUNSPEC(T9) + !ISUNSPEC(T10);

		ADDR(0, 0);
		#if ARGINFO
			typedef arginfo T99[nargs];
			ADDR(end0, 99);
			ADDR(end99, 1);
		#else
			ADDR(end0, 1);
		#endif
		ADDR(end1, 2);
		ADDR(end2, 3);
		ADDR(end3, 4);
		ADDR(end4, 5);
		ADDR(end5, 6);
		ADDR(end6, 7);
		ADDR(end7, 8);
		ADDR(end8, 9);
		ADDR(end9, 10);

		static const int  end = end10;

	protected:
		static const int len = end;				/* length of the packet, assuming last specified variable was a scalar */
		static const int lenp = ASTART(end, ALIGNOF(Tmaxalign));		/* padded length, assuming last variable was a scalar, that properly aligns the next packet */

	public:
		template<typename T>
		__host__ __device__ inline static int len_with_padding(const T& x)		/* padded length of the record when the last element is not an array<> */
		{
			return lenp;
		}

		template<typename T>
		__host__ __device__ inline static int len_with_padding(const array<T> &x)	/* padded length of the record when the last element is an array<> */
		{
			// compute the end offset
			int at2 = end;
			at2 += (x.nelem-1) * SIZEOF(T);

			// add padding to next packet
			int lenp = ASTART(at2, ALIGNOF(Tmaxalign));

			return lenp;
		}

	public:
		#if ARGINFO
		template<typename T>
		__host__ __device__ inline static void store_arginfo(const char *ptr, const T &x)
		{
			return store_arginfo_aux(ptr, 1);
		}

		template<typename T>
		__host__ __device__ inline static void store_arginfo(const char *ptr, const array<T> &x)
		{
			return store_arginfo_aux(ptr, x.nelem);
		}

		__host__ __device__ inline static void store_arginfo_aux(const char *ptr, int nelem)
		{
			DHOST( std::cerr << "Storing arginfo [nargs=" << nargs << " nelem=" << nelem << "]\n"; )
			DHOST( std::cerr << "Unspecified args: " << ISUNSPEC(T0) << " " << ISUNSPEC(T1) << " " << ISUNSPEC(T2) << " " << ISUNSPEC(T3) << " " << ISUNSPEC(T4) << " " << ISUNSPEC(T5) << " " << ISUNSPEC(T6) << " " << ISUNSPEC(T7) << " " << ISUNSPEC(T8) << " " << ISUNSPEC(T9) << " " << ISUNSPEC(T10) << "\n"; )

			arginfo *ai = (arginfo *)(ptr + begin99);
			get_arginfo0(ai + 0, nelem);
			get_arginfo1(ai + 1, nelem);
			get_arginfo2(ai + 2, nelem);
			get_arginfo3(ai + 3, nelem);
			get_arginfo4(ai + 4, nelem);
			get_arginfo5(ai + 5, nelem);
			get_arginfo6(ai + 6, nelem);
			get_arginfo7(ai + 7, nelem);
			get_arginfo8(ai + 8, nelem);
			get_arginfo9(ai + 9, nelem);
			get_arginfo10(ai + 10, nelem);

			// store the data into the header (assumed to be T0)
			T0 &hdr = *(T0 *)(ptr + begin0);
			hdr.nargs = nargs;
			hdr.infos = begin99;
		}
		#endif

		__host__ __device__ inline static void dump()	/* debugging */
		{
			dump0();
		#if ARGINFO
			dump99();
		#endif
			dump1(); dump2(); dump3(); dump4(); dump5(); dump6(); dump7(); dump8(); dump9(); dump10();
		}
	};


	}
}

#endif // bits_gpulog_msg_layout_h__
