#ifndef swarmlog_h__
#define swarmlog_h__

#include "gpulog/gpulog.h"
#include "gpulog/lprintf.h"

#include "swarm.h"

#if __CUDACC__
// The assumption is all CUDA code will be concatenated/included and compiled
// as a single source file (thus avoiding the creation of duplicate copies of 
// hlog and dlog)
gpulog::host_log hlog;
__constant__ gpulog::device_log dlog;
#endif

// declaration for g++-compiled sources
extern gpulog::host_log hlog;

namespace swarm
{
	void init_logs(const std::string &writer_cfg);
	void flush_logs(bool ifneeded = false);
	void sort_binary_log_file(const std::string &outfn, const std::string &infn);

	static const int EVT_SNAPSHOT		= 1000000;	// marks a snapshot of a system. see system_snapshot() down below

	template<typename L, typename T1>
		__host__ __device__ inline PTR_T(SCALAR(T1)) log_event(L &l, const int recid, const double T, const int sys, const T1 &v1)
		{
			return l.write(recid, T, sys, v1);
		}

	template<typename L, typename T1, typename T2>
		__host__ __device__ inline PTR_T(SCALAR(T2)) log_event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2)
		{
			return l.write(recid, T, sys, v1, v2);
		}

	template<typename L, typename T1, typename T2, typename T3>
		__host__ __device__ inline PTR_T(SCALAR(T3)) log_event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3)
		{
			return l.write(recid, T, sys, v1, v2, v3);
		}

	template<typename L, typename T1, typename T2, typename T3, typename T4>
		__host__ __device__ inline PTR_T(SCALAR(T4)) log_event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4)
		{
			return l.write(recid, T, sys, v1, v2, v3, v4);
		}

	template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5>
		__host__ __device__ inline PTR_T(SCALAR(T5)) log_event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5)
		{
			return l.write(recid, T, sys, v1, v2, v3, v4, v5);
		}

	template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
		__host__ __device__ inline PTR_T(SCALAR(T6)) log_event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6)
		{
			return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6);
		}

	template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
		__host__ __device__ inline PTR_T(SCALAR(T7)) log_event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7)
		{
			return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6, v7);
		}

	template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
		__host__ __device__ inline PTR_T(SCALAR(T8)) log_event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8)
		{
			return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6, v7, v8);
		}
	#if 0
	template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
		__host__ __device__ inline PTR_T(SCALAR(T9)) log_event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8, const T9 &v9)
		{
			return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6, v7, v8, v9);
		}

	template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
		__host__ __device__ inline PTR_T(SCALAR(T10)) log_event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8, const T9 &v9, const T10 &v10)
		{
			return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10);
		}
	#endif

	// for on-GPU state logging of bodies
	// NOTE: I've written out the datatypes _explicitly_, because
	// of alignment requirements that have to be hand-tuned between
	// the device and host code. Yes, this _is_ unfortunate.
	struct ALIGN(8) body
	{
		// NOTE: put all doubles first, to avoid interstitial padding
		// and alignment nvcc vs. gcc issues
		double	x, y, z;
		double	vx, vy, vz;

		float	m;
		int bod;

		// load body information from ensemble to body structure
		__device__ __host__ void set(const ensemble &ens, int sys, int bod_)
		{
			bod = bod_;
			m = ens.mass(sys, bod);
			x = ens.x(sys, bod);
			y = ens.y(sys, bod);
			z = ens.z(sys, bod);
			vx = ens.vx(sys, bod);
			vy = ens.vy(sys, bod);
			vz = ens.vz(sys, bod);
		}
	};

	template<int N>
	struct body_set_cls
	{
		const ensemble &ens;
		int sys, bod[N];

		__device__ __host__ inline body_set_cls(const ensemble &ens_, int sys_) : ens(ens_), sys(sys_) { }
	};
	__device__ __host__ inline const body_set_cls<1> body_set(const ensemble &ens, int sys, int bod0)
	{
		body_set_cls<1> br(ens, sys);
		br.bod[0] = bod0;
		return br;
	}
	__device__ __host__ inline const body_set_cls<2> body_set(const ensemble &ens, int sys, int bod0, int bod1)
	{
		body_set_cls<2> br(ens, sys);
		br.bod[0] = bod0;
		br.bod[1] = bod1;
		return br;
	}
	__device__ __host__ inline const body_set_cls<3> body_set(const ensemble &ens, int sys, int bod0, int bod1, int bod2)
	{
		body_set_cls<3> br(ens, sys);
		br.bod[0] = bod0;
		br.bod[1] = bod1;
		br.bod[2] = bod2;
		return br;
	}

	/*
		Store a snapshot of the entire system (EVT_SNAPSHOT).
	*/
	template<typename L>
	__device__ __host__ void log_snapshot(L &l, const ensemble &ens, const int sys, const double T)
	{
		body *bodies = log_event(l, EVT_SNAPSHOT, T, sys, gpulog::array<body>(ens.nbod()));
		for(int bod=0; bod != ens.nbod(); bod++)
		{
			bodies[bod].set(ens, sys, bod);
		}
	}
}

namespace gpulog
{
	namespace internal
	{
		// body_set_cls is a proxy for an array of bodies, so make sure it reports
		// the same alignment as body[N], as well as sizeof()
		template<int N> struct alignment<swarm::body_set_cls<N> > : public alignment<swarm::body[N]> { };	// return alignment of body[N]
		template<int N> struct    ttrait<swarm::body_set_cls<N> > : public    ttrait<swarm::body[N]> { };	// return traits of body[N]

		// serialization of a list of bodies
		template<int N> struct     argio<swarm::body_set_cls<N> >
		{
			__host__ __device__ static inline void put(char *ptr, const swarm::body_set_cls<N> &br, int start, int datalen)
			{
				DHOST( std::cerr << "Writing [" << br << "] start=" << start << " len=" << datalen << "\n" );
				DGPU( printf("Writing start=%d len=%d\n", start, datalen); );
				dev_assert(sizeof(swarm::body)*N == datalen);

				// write out N bodies
				swarm::body *bodies = (swarm::body *)(ptr + start);
				for(int i=0; i != N; i++)
				{
					bodies[i].set(br.ens, br.sys, br.bod[i]);
				}
			}
		};
	}
}

#if 0
#include "swarm.h"
#include <cux/cux.h>
#include <cuda.h>

/**

	Event logging and output facility for swarm

	Before peeking at the horror below, read the description in docs/eventlog.txt
*/


namespace swarm {

// System-defined event IDs. The user must ensure there's no collision between their event IDs and these.
static const int EVT_EOF	= 0;

static const int EVT_PRINTF		= 1000000;	// printf event (mostly good for debugging)
static const int EVT_MSGLOST		= 1000001;	// marker that a message was dropped due to being too long to fit into MAX_MSG_LEN
static const int EVT_SNAPSHOT		= 1000002;	// marks a snapshot of a system. see output_if_needed() in swarmlib.cu
static const int EVT_SNAPSHOT_MARKER	= 1000003;	// a dummy event used to generate a unique snapshot ID. Should be ignored by readers.

extern "C" void debug_hook();

/*
	Align the pointer to a 4-byte boundary, preparing it to store
	the size of next data subpacket.
*/
inline __device__ __host__ void align_for_header(int &at)
{
	while (at % 4) at++;
}

/*
	The function below serves to align the pointer to either
	a 4, 8 or 16 byte boundary, depending on whether the structure
	being written is bigger than 4, 8 or 16 bytes, respectively.
	This is NECESSARY, as the GPU hardware appears to be incapable
	of writing to unaligned locations (e.g., an attempt to write
	a double to a location that is not 8-byte aligned will result
	in silent overwriting of data before/after the location).
*/
inline __device__ __host__ void align_for_payload(int &at, int size)
{
	if(size > 4) { while (at % 8) at++; }
	if(size > 8) { while (at % 16) at++; }
}

// event::SIZEOF long structure, consisting of a header describing the event
// and the body containing the user-supplied event data
struct event // TODO: Rename it event_packet
{
public:
	struct ALIGN(16) header
	{
		int evtid;	// event ID
		int evtref;	// unique event number
		int len;	// payload length
		int threadId;	// originator thread

		header(int evtid_ = 0, int evtref_ = 0, int len_ = 0, int threadId_ = 0)
			: evtid(evtid_), evtref(evtref_), len(len_), threadId(threadId_)
		{ }
	};

	static const int SIZEOF = 256;					// event structure size (in bytes)
	static const int PAYLOAD_SIZE = SIZEOF-sizeof(event::header);	// user-data payload size

	// accessors
	int evtid() const { return hdr.evtid; }		// event ID of the current event
	int evtref() const { return hdr.evtref; }	// unique reference to this event
	int threadId() const { return hdr.threadId; }	// originator threadId of the current event
	int len() const { return hdr.len; }		// byte-length of payload

protected:
	header hdr;
	char	data[PAYLOAD_SIZE];

	friend class ievent;
	friend class oevent;
	friend class eventlog;
};

//
// event data unserializer
//
struct ievent
{
public:
	struct opaque_data;

protected:
	int at;			// offset within evt->data
	const event *evt;

	bool test_end();
public:
	// accessors
	int evtid() const { return evt->hdr.evtid; }		// event ID of the current event
	int evtref() const { return evt->hdr.evtref; }		// unique reference to this event
	int threadId() const { return evt->hdr.threadId; }	// originator threadId of the current event

	ievent(const event &ev) : evt(&ev), at(0) {}
	//ievent() : evt(NULL), at(-1) {}

	ievent &operator >>(opaque_data &v);
	ievent &operator >>(std::string &s);
	ievent &operator >>(char *s);
	ievent &operator >>(event &evt);
	template<typename T> ievent &operator >>(T &v)
	{
		if(!test_end()) { return *this; }

		align_for_header(at);
		int size = *(int*)(evt->data + at);
		if(size != sizeof(T)) ERROR("Programmer error: data size != type size. Contact the authors.");
		at += sizeof(int);

		align_for_payload(at, size);
		v = *(T*)(evt->data + at);
		at += sizeof(T);
		return *this;
	}

	std::string printf() const;
	bool isprintf() const;

	operator bool() const { return evt && at != -1; }
};

class device_eventlog;
class host_eventlog;

//
// event data serializer
//
struct oevent
{
protected:
	int offs;		// offset within evt->data
	event *evt;		// output event packet

public:
	__device__ __host__ int evtid() const { return evt->hdr.evtid; }		// event ID of the current event
	__device__ __host__ int evtref() const { return evt->hdr.evtref; }		// unique reference to this event
	__device__ __host__ int threadId() const { return evt->hdr.threadId; }	// originator threadId of the current event

#ifdef __CUDACC__
	__device__ oevent(device_eventlog &dest, int evtid);
#endif
	__host__ oevent(host_eventlog &dest, int evtid);
	__device__ __host__ ~oevent()
	{
		evt->hdr.len = offs;
	}

	template<typename T>
	__device__ __host__ oevent &operator<<(const T &data)
	{
		align_for_header(offs);
		if(sizeof(int)+offs > event::PAYLOAD_SIZE) { offs = -1; return *this; } // buffer overflow; abort writing.
		*(int*)(evt->data + offs) = sizeof(data); offs += sizeof(int);

		align_for_payload(offs, sizeof(T));
		if(sizeof(T)+offs > event::PAYLOAD_SIZE) { offs = -1; return *this; } // buffer overflow; abort writing.
		*(T*)(evt->data + offs) = data; offs += sizeof(T);

		return *this;
	}

	__device__ __host__ oevent &operator<<(const char *c)
	{
		align_for_header(offs);
		int offs0 = offs;
		offs += sizeof(int);

		// copy the string, including '\0' to the output
		do
		{
			if(offs >= event::PAYLOAD_SIZE)
			{
				// buffer overflow; abort writing, signal overflow.
				offs = -1;
				return *this;
			};
	
			evt->data[offs] = *c;
			offs++;
		} while(*(c++) != '\0');

		// store string length as negative number (signaling it's unaligned)
		*(int*)(evt->data + offs0) = -(offs - offs0 - sizeof(int));

		return *this;
	}
};


// for on-GPU state logging of bodies
// NOTE: I've written out the datatypes _explicitly_, because
// of alignment requirements that have to be hand-tuned between
// the device and host code. Yes, this _is_ unfortunate.
struct ALIGN(8) body
{
	// NOTE: put all doubles first, to avoid interstitial padding
	double	T;
	double	x, y, z;
	double	vx, vy, vz;

	float	m;

	int		sys, bod;
	int		evtref, user_data;	// e.g.: flags, eventId, whatever...
};

// argument-passing struct for writer::process()
struct output_buffers
{
	const event *events;
	int nevt, nevt_dropped;

	const body *bodies;
	int nbod, nbod_dropped;
};

//
// eventlog - Base class for event log objects. DO NOT USE this object directly;
// use instead the derived host_/device_eventlog classes, and hlog/dlog static
// objects.
//
struct eventlog
{
protected:
	//
	// NOTE: This method should _exclusively_ be called by oevent() constructor,
	//       and noone else.
	//
	__device__ __host__ event *prepare_event(int idx, int evtid, int threadid)
	{
		if(idx >= ecap) { return NULL; }
		event *evt = &events.evt[idx];

		evt->hdr.evtref = evtref_base + idx;
		evt->hdr.evtid = evtid;
		evt->hdr.threadId = threadid;
		// NOTE: len MUST be filled out by the caller.

		return evt;
	}

	__device__ __host__ int store_body(int nbod, const ensemble &ens, int sys, int bod, double T, int evtref, int user_data = -1)
	{
		if(nbod >= bcap) { return -1; }
	
		body &b = bodies[nbod];
		b.sys = sys;
		b.bod = bod;
		b.T = T;
		b.m = ens.mass(sys, bod);
		b.x = ens.x(sys, bod);
		b.y = ens.y(sys, bod);
		b.z = ens.z(sys, bod);
		b.vx = ens.vx(sys, bod);
		b.vy = ens.vy(sys, bod);
		b.vz = ens.vz(sys, bod);
		b.evtref = evtref;
		b.user_data = user_data;
	
		return bodref_base + nbod;
	}

public:

	struct counters
	{
		int nbod, nevt;	// current buffer positions

		__device__ __host__ void reset() { nbod = nevt = 0; }
	};

	// buffer capacities and tresholds to trigger flushing
	int bcap, ecap;
	int btresh, etresh;

	// ref bases
	int evtref_base;
	int bodref_base;

	// data
	counters *ctr;
	union
	{
		char *raw;
		event *evt;
	} events;

	//
	// Data/functions specific to integrator code
	//
	body *bodies;
	//system *systems;
};

#if __CUDACC__
struct device_eventlog : public eventlog
{
public:
	__device__ event *alloc_event(int evtid)
	{
		int idx = atomicAdd(&ctr->nevt, 1);
		return prepare_event(idx, evtid, threadId());
	}

	__device__ int log_body(const ensemble &ens, int sys, int bod, double T, int user_data = -1)
	{
		int nbod = atomicAdd(&this->ctr->nbod, 1);
		return store_body(nbod, ens, sys, bod, T, user_data);
	}

public:
	#define __DEVICE__ __device__
	#include "swarmlog_impl.h"
	#undef __DEVICE__
};
#endif

// CPU end of the event log -- receives the events recorded
// on the GPU
struct host_eventlog : public eventlog
{
protected:
	friend class oevent;
	event *alloc_event(int evtid);
	int log_body(const ensemble &ens, int sys, int bod, double T, int user_data = -1);

public:
	#define __DEVICE__
	#include "swarmlog_impl.h"
	#undef __DEVICE__

protected:
	eventlog dlog;	// object used to upload/download the GPU instance
	writer *w;	// sink for events

	bool need_gpu_flush();
	void copyFromGPU();

	bool lastongpu;		// was the last buffer write done on GPU (i.e., have you called prepare_for_gpu())

	void flush_to_writer();

	void prepare_for_cpu();
public:
	void prepare_for_gpu();	// prepare for execution on the GPU. MUST call this function before launching a kernel.

	// initialize GPU and CPU buffers (capacities)
	void initialize(int ecap = 16*1024, int bcap = 1024*1024, int scap = 1024*1024);

	// attach the sink
	void attach_sink(writer *w_);

	// flush the data to sink if CPU/GPU buffers are close to capacity
	void flush_if_needed(bool cpuonly = false);

	// flush the data to sink, syncing with GPU if needed
	void flush();

	// provide access to events buffer
	const event *get_events_buffer(int &nevt) const;
	// provide access to bodies buffer
	const body *get_body_buffer(int &nevt) const;
public:
	host_eventlog();
	~host_eventlog();
};
extern host_eventlog hlog;

#ifdef __CUDACC__
inline __device__ oevent::oevent(device_eventlog &dest, int evtid)
{
	evt = dest.alloc_event(evtid);
	offs = 0;
}
#endif
inline __host__ oevent::oevent(host_eventlog &dest, int evtid)
{
	evt = dest.alloc_event(evtid);
	offs = 0;
}

} // end namespace swarm

#if __CUDACC__
__constant__ swarm::device_eventlog dlog;
#endif

#endif

#endif
