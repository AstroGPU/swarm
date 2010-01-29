#ifndef swarmlog_h__
#define swarmlog_h__

#include "swarm.h"
#include <cux/cux.h>
#include <cuda.h>

/**

	Event logging and output facility for swarm

	Before peeking at the horror below, read the description in docs/eventlog.txt
*/


// System-defined event IDs. The user must ensure there's no collision between their event IDs and these.
static const int EVT_EOF	= 0;

static const int EVT_PRINTF	= 1000000;	// printf event (mostly good for debugging)
static const int EVT_MSGLOST	= 1000001;	// marker that a message was dropped due to being too long to fit into MAX_MSG_LEN
static const int EVT_SNAPSHOT	= 1000002;	// marks a snapshot of a system. see output_if_needed() in swarmlib.cu

struct eventlog_base
{
public:
	static const int MAX_MSG_LEN = 256;	// maximum message length (in bytes)

	struct ALIGN(4) evt_hdr
	{
		int evtid;	// event ID
		int evtref;	// unique event number
		int len;	// payload length
		int threadId;	// originator thread
		int nargs;	// number of sub-packets in the payload

		evt_hdr(int evtid_ = 0, int evtref_ = 0, int len_ = 0, int threadId_ = 0, int nargs_ = 0)
			: evtid(evtid_), evtref(evtref_), len(len_), threadId(threadId_), nargs(nargs_)
		{ }
	};

public:
	struct event
	{
		evt_hdr hdr;
		char	data[MAX_MSG_LEN-sizeof(evt_hdr)];
	};

	struct body // for on-GPU state logging
	{
		int		sys, bod;
		real_time	T;
		real_mass	m;
		real_pos	x, y, z;
		real_vel	vx, vy, vz;
		int		user_data;	// e.g.: flags, eventId, etc.
	};

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
	char *events;

	//
	// Data/functions specific to integrator code
	//
	body *bodies;
	//system *systems;
};

#if __CUDACC__
struct eventlog : public eventlog_base
{
public:
	__device__ void prepare() { }

public:
	#define __DEVICE__ __device__
	#include "swarmlog_impl.h"
	#undef __DEVICE__

	friend struct cpueventlog;
	friend struct ieventstream;
};
__constant__ eventlog glog;
#endif

// CPU end of the event log -- receives the events recorded
// on the GPU
struct cpu_eventlog : public eventlog_base
{
protected:
	int atomicAdd(int *p, int v) { int tmp = *p; *p += v; return tmp; }
	int threadId() { return -1; }
	void prepare()
	{
		prepare_for_cpu();
		flush_if_needed(true);
	}

public:
	#define __DEVICE__
	#include "swarmlog_impl.h"
	#undef __DEVICE__

protected:
	eventlog_base glog;	// object used to upload/download the GPU instance
	writer *w;		// sink for events

	bool need_gpu_flush();
	void copyFromGPU();

	bool lastongpu;		// was the last buffer write done on GPU (i.e., have you called prepare_for_gpu())

	void prepare_for_cpu();
public:
	void prepare_for_gpu();	// prepare for execution on the GPU. MUST call this function before launching a kernel.

	// initialize GPU and CPU buffers (capacities)
	void initialize(int ecap = 16*1024, int bcap = 1024*1024, int scap = 1024*1024);

	// attach the sink
	void attach_sink(writer *w_) { w = w_; }

	// flush the data to sink if CPU/GPU buffers are close to capacity
	void flush_if_needed(bool cpuonly = false);

	// flush the data to sink, syncing with GPU if needed
	void flush();

public:
	cpu_eventlog();
	~cpu_eventlog();
};
extern cpu_eventlog clog;

void debug_hook();

// event stream -- stream-like interface to extract events from cpu_eventlog
struct ieventstream
{
protected:
	friend bool get_message(ieventstream &evt, std::string &res);

protected:
	struct raw_evt
	{
		int size;
		const char *data;
	};

	typedef eventlog_base::evt_hdr evt_hdr;
	typedef eventlog_base::body body;
	typedef eventlog_base::event event;

	evt_hdr hdr;
	const char *data;
	int nevt, m_ndropped;
	int at, atevt;
	bool m_eom;		// end-of-message flag

	cpu_eventlog &log;	// the log to which we're bound

public:
	int evtid() const { return hdr.evtid; }			// event ID of the current event
	int evtref() const { return hdr.evtref; }		// unique reference to this event
	int threadId() const { return hdr.threadId; }		// originator threadId of the current event
	int nargs() const { return hdr.nargs; }			// number of data elements in the current event

	int nevents_dropped() const { return m_ndropped; }		// number of events that were dropped due to GPU buffer exhaustion

	const body* get_bodies(int &nbod, int &ndropped) const;	// get a pointer to recorded bodies

	// advance to next message, returning its event id
	// if called immediately after construction, it will move to first message
	int next();

	void init(const char *data_, int nevt_, int ecap_);
	ieventstream(cpu_eventlog &log);

	bool check_end();
	ieventstream &operator >>(raw_evt &v);
	ieventstream &operator >>(std::string &s);
	ieventstream &operator >>(char *s);
	ieventstream &operator >>(event &evt);
	template<typename T> ieventstream &operator >>(T &v)
	{
		if(!check_end()) { return *this; }

		int size = *(int*)(data + at);
		if(size != sizeof(T)) ERROR("Programmer error: data size != type size. Contact the authors.");
		at += sizeof(int);
		v = *(T*)(data+at);
		at += sizeof(T);
		return *this;
	}

	operator bool() const { return !eof() && !m_eom; }
	bool eof() const { return atevt == nevt; }

	// find and return the next printf message
	std::string next_message();
};

#endif
