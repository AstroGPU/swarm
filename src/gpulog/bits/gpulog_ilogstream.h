#ifndef bits_gpulog_ilogstream_h__
#define bits_gpulog_ilogstream_h__

namespace gpulog
{
namespace internal
{

	/* a stream of logrecords */
	struct ilogstream
	{
	protected:
		header hend;

	protected:
		const char *const ptr;
		int at, len;

	public:
		__device__ __host__ ilogstream(const char *ptr_, int len_) : ptr(ptr_), at(0), len(len_), hend(-1, 0) {}
		__device__ __host__ ilogstream(const host_log &l) : ptr(l.internal_buffer()), at(0), len(l.size()), hend(-1, 0) {}

		__device__ __host__ logrecord next()
		{
			// if we're at the end
			if(at == len) { return logrecord((char *)&hend); }

			// return next packet
			logrecord rec(ptr + at);
			at += rec.len();

			return rec;
		}
	};

}
}

#endif
