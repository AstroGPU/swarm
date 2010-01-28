#ifndef swarmlog_h__
#define swarmlog_h__

#include "swarm.h"
#include <cux/cux.h>
#include <cuda.h>

struct eventlog
{
public:
	static const int MAX_MSG_LEN = 256;	// maximum message length (in bytes)

	struct ALIGN(4) evt_hdr
	{
		int evtid;	// event ID
		int len;	// payload length
		int threadId;	// originator thread
		int nargs;	// number of sub-packets in the payload

		evt_hdr(int evtid_ = 0, int len_ = 0, int threadId_ = 0, int nargs_ = 0)
			: evtid(evtid_), len(len_), threadId(threadId_), nargs(nargs_)
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

	static const int EVT_EOF	= 0;
	static const int EVT_PRINTF	= 999999;
	static const int EVT_MSGLOST	= 999998;

protected:
	struct counters
	{
		int nbodX, nevtX;	// current buffer position
		int needflush;		// a kernel will set this flag to request buffer flushing

		__device__ __host__ void reset() { nbodX = nevtX = 0; needflush = 0; }
	};

	// buffer capacities
	int bcapX, ecapX;
	int btresh, etresh;

	// data
	counters *ctr;
	char *events;

public:
	int needflush() const { return ctr->needflush; }

	//
	// General event logging facility
	//
	template<typename T1>
		__device__ void log_event(int evtId, const T1 &v1);
	template<typename T1, typename T2>
		__device__ void log_event(int evtId, const T1 &v1, const T2 &v2);
	template<typename T1, typename T2, typename T3>
		__device__ void log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3);
	template<typename T1, typename T2, typename T3, typename T4>
		__device__ void log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4);
	template<typename T1, typename T2, typename T3, typename T4, typename T5>
		__device__ void log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5);
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
		__device__ void log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6);
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
		__device__ void log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7);
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
		__device__ void log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8);
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
		__device__ void log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8, const T8 &v8);

	//
	// Printf-like facility
	//
	template<typename T1>
		__device__ void printf(const char *fmt, const T1 &v1)
		{ log_event(EVT_PRINTF, fmt, v1); }
	template<typename T1, typename T2>
		__device__ void printf(const char *fmt, const T1 &v1, const T2 &v2)
		{ log_event(EVT_PRINTF, fmt, v1, v2); }
	template<typename T1, typename T2, typename T3>
		__device__ void printf(const char *fmt, const T1 &v1, const T2 &v2, const T3 &v3)
		{ log_event(EVT_PRINTF, fmt, v1, v2, v3); }
	template<typename T1, typename T2, typename T3, typename T4>
		__device__ void printf(const char *fmt, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4)
		{ log_event(EVT_PRINTF, fmt, v1, v2, v3, v4); }
	template<typename T1, typename T2, typename T3, typename T4, typename T5>
		__device__ void printf(const char *fmt, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5)
		{ log_event(EVT_PRINTF, fmt, v1, v2, v3, v4, v5); }
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
		__device__ void printf(const char *fmt, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6)
		{ log_event(EVT_PRINTF, fmt, v1, v2, v3, v4, v5, v6); }
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
		__device__ void printf(const char *fmt, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7)
		{ log_event(EVT_PRINTF, fmt, v1, v2, v3, v4, v5, v6, v7); }
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
		__device__ void printf(const char *fmt, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8)
		{ log_event(EVT_PRINTF, fmt, v1, v2, v3, v4, v5, v6, v7, v8); }
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
		__device__ void printf(const char *fmt, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8, const T9 &v9)
		{ log_event(EVT_PRINTF, fmt, v1, v2, v3, v4, v5, v6, v7, v8, v9); }

protected:
	template<typename T>
		__device__ void push_data(int &nevt, int mend, const T &data);
	__device__ void push_data(int &nevt, int mend, const char *c);
	__device__ bool evt_start(int &nevt, int &mend);
	__device__ void evt_end(int &nevt, int &mend, int nargs, int evtid);

	friend struct cpueventlog;
	friend struct ieventstream;

	//
	// Data/functions specific to integrator code
	//
protected:
	body *bodies;
	//system *systems;

public:
	__device__ bool log_system(const ensemble &ens, int sys, double T, int user_data = -1);
	__device__ bool log_body(const ensemble &ens, int sys, int bod, double T, int user_data = -1);
};

// event stream
struct ieventstream
{
protected:
	friend bool get_as_message(ieventstream &msg, std::string &res);

protected:
	struct raw_msg
	{
		int size;
		const char *data;
	};

	eventlog::evt_hdr hdr;
	const char *data;
	int nmsg, m_nlost;
	int at, atmsg;
	bool m_eom;		// end-of-message flag

public:
	int evtid() const { return hdr.evtid; }			// event ID of the current event
	int threadId() const { return hdr.threadId; }		// originator threadId of the current event
	int nargs() const { return hdr.nargs; }			// number of data elements in the current event

	int nevents_dropped() const { return m_nlost; }		// number of events that were dropped due to GPU buffer exhaustion

	const eventlog::body* get_bodies(int &nbod, int &ndropped) const;	// get a pointer to recorded bodies

	// advance to next message, returning its event id
	// if called immediately after construction, it will move to first message
	int next();

	void init(const char *data_, int nmsg_, int ecap_);
	ieventstream(bool syncNew = true);

	bool check_end();
	ieventstream &operator >>(raw_msg &v);
	ieventstream &operator >>(std::string &s);
	ieventstream &operator >>(char *s);
	ieventstream &operator >>(eventlog::event &evt);
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
	bool eof() const { return atmsg == nmsg; }

	// find and return the next printf message
	std::string next_message();
};

void debug_hook();
void initialize_eventlog(int ecap = 64, int bcap = 1024*1024, int scap = 1024*1024);

#if __CUDACC__

__constant__ eventlog glog;

__device__ bool eventlog::log_body(const ensemble &ens, int sys, int bod, double T, int user_data)
{
	int nbod = atomicAdd(&ctr->nbodX, 1);
	if(nbod >= bcapX) { return false; }
	if(nbod == btresh) { ctr->needflush = 1; }

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
	b.user_data = user_data;

	return true;
}

__device__ bool eventlog::log_system(const ensemble &ens, int sys, double T, int user_data)
{
	bool succ = true;
	for(int bod = 0; bod != ens.nbod(); bod++)
	{
		succ &= log_body(ens, sys, bod, T);
	}
	return succ;
}

template<typename T>
__device__ void eventlog::push_data(int &nevt, int mend, const T &data)
{
	if(sizeof(int)+sizeof(T)+nevt > mend)
	{
		// buffer overflow; abort writing.
		nevt = mend+1;
		return;
	}

	*(int*)(events + nevt) = sizeof(data); nevt += sizeof(int);
	*(T*)(events + nevt) = data; nevt += sizeof(T);
}

__device__ void eventlog::push_data(int &nevt, int mend, const char *c)
{
	int nevt0 = nevt;
	nevt += sizeof(int);

	// copy the string, including '\0' to the output
	do
	{
		if(nevt >= mend)
		{
			// buffer overflow; abort writing, signal overflow.
			nevt = mend+1;
			return;
		};

		events[nevt] = *c;
		nevt++;
	} while(*(c++) != '\0');

	// store string length
	*(int*)(events + nevt0) = nevt - nevt0 - sizeof(int);
}

__device__ bool eventlog::evt_start(int &nevt, int &mend)
{
	nevt = atomicAdd(&ctr->nevtX, 1);
	if(nevt >= ecapX) { return false; }
	
	// request flushing if we're getting close to the
	// end of the buffer
	if(nevt == etresh)
	{
		ctr->needflush = 1;
	}

	nevt *= MAX_MSG_LEN;
	mend  = nevt + MAX_MSG_LEN;

	// leave room for evt_hdr, which will be written
	// by evt_end
	nevt += sizeof(evt_hdr);

	return true;
}

__device__ void eventlog::evt_end(int &nevt, int &mend, int nargs, int evtid)
{
	if(nevt > mend)
	{
		// rewind and store which message was lost
		nevt = mend - MAX_MSG_LEN + sizeof(evt_hdr);
		push_data(nevt, mend, evtid);
		evtid = EVT_MSGLOST;
	}

	// store the event header
	evt_hdr hdr(evtid, nevt-(mend-MAX_MSG_LEN), threadId(), nargs);
	nevt = mend - MAX_MSG_LEN;
	*(evt_hdr*)(events + nevt) = hdr;
}

#define MSG_START(evtid, nargs) \
	int mend, nevt; \
	if(!evt_start(nevt, mend)) { return; }

#define MSG_END(evtid, nargs) \
	evt_end(nevt, mend, nargs, evtid)

template<typename T1>
	__device__ void eventlog::log_event(int evtId, const T1 &v1)
{
	MSG_START(evtId, 1);

	push_data(nevt, mend, v1);

	MSG_END(evtId, 1);
}


template<typename T1, typename T2>
	__device__ void eventlog::log_event(int evtId, const T1 &v1, const T2 &v2)
{
	MSG_START(evtId, 2);

	push_data(nevt, mend, v1);
	push_data(nevt, mend, v2);

	MSG_END(evtId, 2);
}
template<typename T1, typename T2, typename T3>
	__device__ void eventlog::log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3)
{
	MSG_START(evtId, 3);

	push_data(nevt, mend, v1);
	push_data(nevt, mend, v2);
	push_data(nevt, mend, v3);

	MSG_END(evtId, 3);
}
template<typename T1, typename T2, typename T3, typename T4>
	__device__ void eventlog::log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4)
{
	MSG_START(evtId, 4);

	push_data(nevt, mend, v1);
	push_data(nevt, mend, v2);
	push_data(nevt, mend, v3);
	push_data(nevt, mend, v4);

	MSG_END(evtId, 4);
}
template<typename T1, typename T2, typename T3, typename T4, typename T5>
	__device__ void eventlog::log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5)
{
	MSG_START(evtId, 5);

	push_data(nevt, mend, v1);
	push_data(nevt, mend, v2);
	push_data(nevt, mend, v3);
	push_data(nevt, mend, v4);
	push_data(nevt, mend, v5);

	MSG_END(evtId, 5);
}
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
	__device__ void eventlog::log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6)
{
	MSG_START(evtId, 6);

	push_data(nevt, mend, v1);
	push_data(nevt, mend, v2);
	push_data(nevt, mend, v3);
	push_data(nevt, mend, v4);
	push_data(nevt, mend, v5);
	push_data(nevt, mend, v6);

	MSG_END(evtId, 6);
}
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
	__device__ void eventlog::log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7)
{
	MSG_START(evtId, 7);

	push_data(nevt, mend, v1);
	push_data(nevt, mend, v2);
	push_data(nevt, mend, v3);
	push_data(nevt, mend, v4);
	push_data(nevt, mend, v5);
	push_data(nevt, mend, v6);
	push_data(nevt, mend, v7);

	MSG_END(evtId, 7);
}
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
	__device__ void eventlog::log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8)
{
	MSG_START(evtId, 8);

	push_data(nevt, mend, v1);
	push_data(nevt, mend, v2);
	push_data(nevt, mend, v3);
	push_data(nevt, mend, v4);
	push_data(nevt, mend, v5);
	push_data(nevt, mend, v6);
	push_data(nevt, mend, v7);
	push_data(nevt, mend, v8);

	MSG_END(evtId, 8);
}
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
	__device__ void eventlog::log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8, const T8 &v8)
{
	MSG_START(evtId, 9);

	push_data(nevt, mend, v1);
	push_data(nevt, mend, v2);
	push_data(nevt, mend, v3);
	push_data(nevt, mend, v4);
	push_data(nevt, mend, v5);
	push_data(nevt, mend, v6);
	push_data(nevt, mend, v7);
	push_data(nevt, mend, v8);
	push_data(nevt, mend, v9);

	MSG_END(evtId, 9);
}

#endif

#endif
