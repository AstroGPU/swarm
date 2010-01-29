	//
	// Printf-like facility
	//
		__DEVICE__ int printf(const char *fmt)
		{ return log_event(EVT_PRINTF, fmt); }
	template<typename T1>
		__DEVICE__ int printf(const char *fmt, const T1 &v1)
		{ return log_event(EVT_PRINTF, fmt, v1); }
	template<typename T1, typename T2>
		__DEVICE__ int printf(const char *fmt, const T1 &v1, const T2 &v2)
		{ return log_event(EVT_PRINTF, fmt, v1, v2); }
	template<typename T1, typename T2, typename T3>
		__DEVICE__ int printf(const char *fmt, const T1 &v1, const T2 &v2, const T3 &v3)
		{ return log_event(EVT_PRINTF, fmt, v1, v2, v3); }
	template<typename T1, typename T2, typename T3, typename T4>
		__DEVICE__ int printf(const char *fmt, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4)
		{ return log_event(EVT_PRINTF, fmt, v1, v2, v3, v4); }
	template<typename T1, typename T2, typename T3, typename T4, typename T5>
		__DEVICE__ int printf(const char *fmt, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5)
		{ return log_event(EVT_PRINTF, fmt, v1, v2, v3, v4, v5); }
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
		__DEVICE__ int printf(const char *fmt, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6)
		{ return log_event(EVT_PRINTF, fmt, v1, v2, v3, v4, v5, v6); }
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
		__DEVICE__ int printf(const char *fmt, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7)
		{ return log_event(EVT_PRINTF, fmt, v1, v2, v3, v4, v5, v6, v7); }
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
		__DEVICE__ int printf(const char *fmt, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8)
		{ return log_event(EVT_PRINTF, fmt, v1, v2, v3, v4, v5, v6, v7, v8); }
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
		__DEVICE__ int printf(const char *fmt, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8, const T9 &v9)
		{ return log_event(EVT_PRINTF, fmt, v1, v2, v3, v4, v5, v6, v7, v8, v9); }

public:
	__DEVICE__ int log_body(const ensemble &ens, int sys, int bod, double T, int user_data = -1)
	{
		this->prepare();
		int nbod = atomicAdd(&this->ctr->nbod, 1);
		if(nbod >= this->bcap) { return nbod + this->bodref_base; }
	
		body &b = this->bodies[nbod];
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
	
		return this->bodref_base + nbod;
	}
	
	__DEVICE__ void log_system(const ensemble &ens, int sys, double T, int user_data = -1)
	{
		for(int bod = 0; bod != ens.nbod(); bod++)
		{
			log_body(ens, sys, bod, T, user_data);
		}
	}
	
	template<typename T>
	__DEVICE__ void push_data(int &nevt, int mend, const T &data)
	{
		align_to_header(nevt);
		if(sizeof(int)+nevt > mend) { nevt = mend+1; return; } // buffer overflow; abort writing.
		*(int*)(this->events + nevt) = sizeof(data); nevt += sizeof(int);

		align_to_payload(nevt, sizeof(T));
		if(sizeof(T)+nevt > mend) { nevt = mend+1; return; } // buffer overflow; abort writing.
		*(T*)(this->events + nevt) = data; nevt += sizeof(T);
	}

	__DEVICE__ void push_data(int &nevt, int mend, const char *c)
	{
		align_to_header(nevt);
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
	
			this->events[nevt] = *c;
			nevt++;
		} while(*(c++) != '\0');

		// store string length as negative number (signaling it's unaligned)
		*(int*)(this->events + nevt0) = -(nevt - nevt0 - sizeof(int));
	}

	__DEVICE__ bool evt_start(int &nevt, int &mend)
	{
		this->prepare();
		nevt = atomicAdd(&this->ctr->nevt, 1);
		if(nevt >= this->ecap) { return false; }
	
		nevt *= eventlog_base::MAX_MSG_LEN;		// convert nevt to byte offset
		mend  = nevt + eventlog_base::MAX_MSG_LEN;
	
		// leave room for evt_hdr, which will be written
		// by evt_end
		nevt += sizeof(eventlog_base::evt_hdr);
	
		return true;
	}
	
	__DEVICE__ int evt_end(int &nevt, int &mend, int nargs, int evtid)
	{
		if(nevt > mend)
		{
			// rewind and store which message was lost
			nevt = mend - eventlog_base::MAX_MSG_LEN + sizeof(eventlog_base::evt_hdr);
			push_data(nevt, mend, evtid);
			evtid = EVT_MSGLOST;
		}

		// store the event header
		int evtref = this->evtref_base;
		evtref += nevt / eventlog_base::MAX_MSG_LEN;
		eventlog_base::evt_hdr hdr(evtid, evtref, nevt-(mend-eventlog_base::MAX_MSG_LEN), threadId(), nargs);
		nevt = mend - eventlog_base::MAX_MSG_LEN;
		*(eventlog_base::evt_hdr*)(this->events + nevt) = hdr;

		return evtref;
	}
	
	#define MSG_START(evtid, nargs) \
		int mend, nevt; \
		if(!evt_start(nevt, mend)) { return nevt; }
	
	#define MSG_END(evtid, nargs) \
		return evt_end(nevt, mend, nargs, evtid)
	
	template<typename T1>
		__DEVICE__ int log_event(int evtId, const T1 &v1)
	{
		MSG_START(evtId, 1);
	
		push_data(nevt, mend, v1);
	
		MSG_END(evtId, 1);
	}
	
	
	template<typename T1, typename T2>
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2)
	{
		MSG_START(evtId, 2);
	
		push_data(nevt, mend, v1);
		push_data(nevt, mend, v2);
	
		MSG_END(evtId, 2);
	}
	template<typename T1, typename T2, typename T3>
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3)
	{
		MSG_START(evtId, 3);
	
		push_data(nevt, mend, v1);
		push_data(nevt, mend, v2);
		push_data(nevt, mend, v3);
	
		MSG_END(evtId, 3);
	}
	template<typename T1, typename T2, typename T3, typename T4>
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4)
	{
		MSG_START(evtId, 4);
	
		push_data(nevt, mend, v1);
		push_data(nevt, mend, v2);
		push_data(nevt, mend, v3);
		push_data(nevt, mend, v4);
	
		MSG_END(evtId, 4);
	}
	template<typename T1, typename T2, typename T3, typename T4, typename T5>
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5)
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
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6)
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
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7)
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
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8)
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
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8, const T9 &v9)
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
