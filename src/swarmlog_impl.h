#if 1
	// NOTE: We could rewrite this to be MUCH, MUCH simpler if we could use
	//       variadic templates from C++0x
	//
		__DEVICE__ int log_event(int evtid)
		{ return (oevent(*this, evtid)).evtref(); }

	template<typename T1>
		__DEVICE__ int log_event(int evtid, const T1 &v1)
		{ return (oevent(*this, evtid) << v1).evtref(); }

	template<typename T1, typename T2>
		__DEVICE__ int log_event(int evtid, const T1 &v1, const T2 &v2)
		{ return (oevent(*this, evtid) << v1 << v2).evtref(); }

	template<typename T1, typename T2, typename T3>
		__DEVICE__ int log_event(int evtid, const T1 &v1, const T2 &v2, const T3 &v3)
		{ return (oevent(*this, evtid) << v1 << v2 << v3).evtref(); }

	template<typename T1, typename T2, typename T3, typename T4>
		__DEVICE__ int log_event(int evtid, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4)
		{ return (oevent(*this, evtid) << v1 << v2 << v3 << v4).evtref(); }

	template<typename T1, typename T2, typename T3, typename T4, typename T5>
		__DEVICE__ int log_event(int evtid, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5)
		{ return (oevent(*this, evtid) << v1 << v2 << v3 << v4 << v5).evtref(); }

	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
		__DEVICE__ int log_event(int evtid, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6)
		{ return (oevent(*this, evtid) << v1 << v2 << v3 << v4 << v5 << v6).evtref(); }

	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
		__DEVICE__ int log_event(int evtid, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7)
		{ return (oevent(*this, evtid) << v1 << v2 << v3 << v4 << v5 << v6 << v7).evtref(); }
	
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
		__DEVICE__ int log_event(int evtid, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8)
		{ return (oevent(*this, evtid) << v1 << v2 << v3 << v4 << v5 << v6 << v7 << v8).evtref(); }
	
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
		__DEVICE__ int log_event(int evtid, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8, const T9 &v9)
		{ return (oevent(*this, evtid) << v1 << v2 << v3 << v4 << v5 << v6 << v7 << v8 << v9).evtref(); }
	
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
#endif

public:	
	__DEVICE__ void log_system(const ensemble &ens, int sys, double T)
	{
		// write out event header
		int snapid = log_event(EVT_SNAPSHOT_MARKER);
		debug_hook();

		int bod_begin, bod_last;
		for(int bod = 0; bod != ens.nbod(); bod++)
		{
			bod_last = log_body(ens, sys, bod, T, snapid);
			if(bod == 0) { bod_begin = bod_last; }
		}

		log_event(EVT_SNAPSHOT, (int)sys, (int)ens.nbod(), (double)T, (int)snapid, (int)bod_begin, (int)bod_last);
	}
#if 0
	template<typename T>
	__DEVICE__ void push_data(int &offs, int endoffs, const T &data)
	{
		align_for_header(offs);
		if(sizeof(int)+offs > endoffs) { offs = endoffs+1; return; } // buffer overflow; abort writing.
		*(int*)(this->events.raw + offs) = sizeof(data); offs += sizeof(int);

		align_for_payload(offs, sizeof(T));
		if(sizeof(T)+offs > endoffs) { offs = endoffs+1; return; } // buffer overflow; abort writing.
		*(T*)(this->events.raw + offs) = data; offs += sizeof(T);
	}

	__DEVICE__ void push_data(int &offs, int endoffs, const char *c)
	{
		align_for_header(offs);
		int offs0 = offs;
		offs += sizeof(int);

		// copy the string, including '\0' to the output
		do
		{
			if(offs >= endoffs)
			{
				// buffer overflow; abort writing, signal overflow.
				offs = endoffs+1;
				return;
			};
	
			this->events.raw[offs] = *c;
			offs++;
		} while(*(c++) != '\0');

		// store string length as negative number (signaling it's unaligned)
		*(int*)(this->events.raw + offs0) = -(offs - offs0 - sizeof(int));
	}

	__DEVICE__ bool evt_start(int &offs, int &endoffs)
	{
		this->prepare();
		offs = atomicAdd(&this->ctr->nevt, 1);
		if(offs >= this->ecap) { return false; }

		offs *= event::SIZEOF;		// convert offs to byte offset
		endoffs  = offs + event::SIZEOF;
	
		// leave room for evt_hdr, which will be written
		// by evt_end
		offs += sizeof(event::header);
	
		return true;
	}
	
	__DEVICE__ int evt_end(int &offs, int &endoffs, int nargs, int evtid)
	{
		if(offs > endoffs)
		{
			// rewind and store which message was lost
			offs = endoffs - event::SIZEOF + sizeof(event::header);
			push_data(offs, endoffs, evtid);
			evtid = EVT_MSGLOST;
		}

		// store the event header
		int evtref = this->evtref_base;
		evtref += offs / event::SIZEOF;
		event::header hdr(evtid, evtref, offs-(endoffs-event::SIZEOF), threadId(), nargs);
		offs = endoffs - event::SIZEOF;
// 		*(event::header*)(this->events.raw + offs) = hdr;

		return evtref;
	}
	
	#define MSG_START(evtid, nargs) \
		int offs, endoffs; \
		if(!evt_start(offs, endoffs)) { return offs; }
	
	#define MSG_END(evtid, nargs) \
		return evt_end(offs, endoffs, nargs, evtid)
	
	template<typename T1>
		__DEVICE__ int log_event(int evtId, const T1 &v1)
	{
		MSG_START(evtId, 1);
	
		push_data(offs, endoffs, v1);
	
		MSG_END(evtId, 1);
	}
	
	
	template<typename T1, typename T2>
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2)
	{
		MSG_START(evtId, 2);
	
		push_data(offs, endoffs, v1);
		push_data(offs, endoffs, v2);
	
		MSG_END(evtId, 2);
	}
	template<typename T1, typename T2, typename T3>
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3)
	{
		MSG_START(evtId, 3);
	
		push_data(offs, endoffs, v1);
		push_data(offs, endoffs, v2);
		push_data(offs, endoffs, v3);
	
		MSG_END(evtId, 3);
	}
	template<typename T1, typename T2, typename T3, typename T4>
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4)
	{
		MSG_START(evtId, 4);
	
		push_data(offs, endoffs, v1);
		push_data(offs, endoffs, v2);
		push_data(offs, endoffs, v3);
		push_data(offs, endoffs, v4);
	
		MSG_END(evtId, 4);
	}
	template<typename T1, typename T2, typename T3, typename T4, typename T5>
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5)
	{
		MSG_START(evtId, 5);
	
		push_data(offs, endoffs, v1);
		push_data(offs, endoffs, v2);
		push_data(offs, endoffs, v3);
		push_data(offs, endoffs, v4);
		push_data(offs, endoffs, v5);
	
		MSG_END(evtId, 5);
	}
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6)
	{
		MSG_START(evtId, 6);
	
		push_data(offs, endoffs, v1);
		push_data(offs, endoffs, v2);
		push_data(offs, endoffs, v3);
		push_data(offs, endoffs, v4);
		push_data(offs, endoffs, v5);
		push_data(offs, endoffs, v6);
	
		MSG_END(evtId, 6);
	}
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7)
	{
		MSG_START(evtId, 7);
	
		push_data(offs, endoffs, v1);
		push_data(offs, endoffs, v2);
		push_data(offs, endoffs, v3);
		push_data(offs, endoffs, v4);
		push_data(offs, endoffs, v5);
		push_data(offs, endoffs, v6);
		push_data(offs, endoffs, v7);
	
		MSG_END(evtId, 7);
	}
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8)
	{
		MSG_START(evtId, 8);
	
		push_data(offs, endoffs, v1);
		push_data(offs, endoffs, v2);
		push_data(offs, endoffs, v3);
		push_data(offs, endoffs, v4);
		push_data(offs, endoffs, v5);
		push_data(offs, endoffs, v6);
		push_data(offs, endoffs, v7);
		push_data(offs, endoffs, v8);
	
		MSG_END(evtId, 8);
	}
	template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
		__DEVICE__ int log_event(int evtId, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8, const T9 &v9)
	{
		MSG_START(evtId, 9);
	
		push_data(offs, endoffs, v1);
		push_data(offs, endoffs, v2);
		push_data(offs, endoffs, v3);
		push_data(offs, endoffs, v4);
		push_data(offs, endoffs, v5);
		push_data(offs, endoffs, v6);
		push_data(offs, endoffs, v7);
		push_data(offs, endoffs, v8);
		push_data(offs, endoffs, v9);
	
		MSG_END(evtId, 9);
	}
#endif
