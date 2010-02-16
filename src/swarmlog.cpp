#include "swarmlog.h"
#include <iostream>
#include <astro/binarystream.h>
#include <fstream>
#include <sstream>
#include <memory>

extern "C" void debug_hook()
{
	std::cerr << "";
}

struct ievent::opaque_data
{
	int size;
	const char *data;
};


// Test if we hit the end of message, and set at = -1 if we did. at = -1
// signals an EOF condition (operator bool() will return false).
bool ievent::test_end()
{
	if(evt->hdr.len == at) { at = -1; }
	return *this;
}

ievent &ievent::operator >>(opaque_data &v)
{
	if(!test_end()) { return *this; }

	align_for_header(at);
	v.size = *(int*)(evt->data + at);
	at += sizeof(int);

	if(v.size > 0) { align_for_payload(at, v.size); } else { v.size *= -1; }
	v.data = evt->data+at;
	at += v.size;
	return *this;
}

ievent &ievent::operator >>(std::string &s)
{
	if(!test_end()) { return *this; }

	align_for_header(at);
	int size = *(int*)(evt->data + at);
	assert(size < 0);
	size *= -1;
	at += sizeof(int);

	s = (evt->data+at);
	int strlen = s.size();
	if(strlen != size-1-sizeof(int))
		ERROR("Programmer error: data size != string size. Contact the authors.");
	at += size;
	return *this;
}

ievent &ievent::operator >>(char *s)
{
	if(!test_end()) { return *this; }

	align_for_header(at);
	int size = *(int*)(evt->data + at);
	assert(size < 0);
	size *= -1;
	at += sizeof(int);

	strcpy(s, evt->data+at);
	int strl = strlen(s);
	if(strl != size-1)
		ERROR("Programmer error: data size != string size. Contact the authors.");
	at += size;
	return *this;
}

ievent &ievent::operator >>(event &evt)
{
	if(!test_end()) { return *this; }
	assert(sizeof(event) == event::SIZEOF);
	memcpy(&evt, this->evt->data, sizeof(event));
	at = evt.hdr.len;

	return *this;
}

event *host_eventlog::alloc_event(int evtid)
{
	prepare_for_cpu();
	int idx = ctr->nevt++;
	return prepare_event(idx, evtid, -1);
}

int host_eventlog::log_body(const ensemble &ens, int sys, int bod, double T, int user_data)
{
	prepare_for_cpu();
	flush_if_needed(true);
	int nbod = ctr->nbod++;
	return store_body(nbod, ens, sys, bod, T, user_data);
}

# if 0
const body* ieventstream::get_bodies(int &nbod, int &ndropped) const
{
	nbod = std::min(log.ctr->nbod, log.bcap);
	ndropped = std::max(log.ctr->nbod - log.bcap, 0);
	return log.bodies;
}

ieventstream::ieventstream(host_eventlog &log_)
	: log(log_)
{
	init(log.events, log.ctr->nevt, log.ecap);
}

// advance to next message
int ieventstream::next()
{
	// already at end?
	if(atevt == nevt) { return EVT_EOF; }

	// move and exit if we're at the end
	atevt++;
	if(atevt != 0)
	{
		data += event::SIZEOF;
	}

	// learn the basics about the current message
	if(atevt != nevt)
	{
		hdr = *(eventlog::evt_hdr*)(data);
		at = sizeof(eventlog::evt_hdr);
		m_eom = false;
	}
	else
	{
		hdr.evtid = EVT_EOF;
	}

	// return current message event ID
	return evtid();
}

void ieventstream::init(const char *data_, int nevt_, int ecap_)
{
	data = data_;
	nevt = std::min(nevt_, ecap_);
	m_ndropped = std::max(nevt_ - ecap_, 0);

	// position us _ahead_ of the first message
	// (expecting a call to next() will move us to first message)
	atevt = -1;
	m_eom = true;
}

#endif

host_eventlog hlog;

bool host_eventlog::need_gpu_flush()
{
	// download GPU counters
	eventlog::counters gctr;
	memcpyToHost(&gctr, dlog.ctr);

	return gctr.nevt >= dlog.ecap || gctr.nbod >= dlog.bcap;
}

void host_eventlog::flush_if_needed(bool cpuonly)
{
	bool cpuneedflush = ctr->nevt >= ecap || ctr->nbod >= bcap;
	bool gpuneedflush = !cpuonly && need_gpu_flush();

	if(cpuneedflush || gpuneedflush)
	{
		flush();
	}
}

// download the GPU log data, overwriting the CPU log
// called _ONLY_ from flush()
void host_eventlog::copyFromGPU()
{
	// download GPU log data
	if(ctr == NULL)
		ERROR("Programmer error: host_eventlog must be initialized before use. Contact the authors.");

	// download GPU counters
	memcpyToHost(ctr, dlog.ctr);

	// download the data
	events.evt = hostAlloc(events.evt, ecap);
	memset(events.evt, 0, ecap * sizeof(event)); // debugging!
	memcpyToHost(events.evt, dlog.events.evt, ecap);
	bodies = hostAlloc(bodies, bcap);
	memcpyToHost(bodies, dlog.bodies, bcap);
}

void host_eventlog::prepare_for_gpu()
{
	if(lastongpu) { return; }

	// last write was on the CPU, and may have moved the
	// refbases. Upload them to the GPU. Because of how
	// we do things, the GPU buffers should be empty at
	// this point (it's a bug if they're not)

#if __DEVICE_EMULATION__
	assert(dlog.ctr->nbod == 0 && dlog.ctr->nevt == 0);
#endif

	// update the GPU copy with the current ref bases
	dlog.evtref_base = evtref_base + ctr->nevt;
	dlog.bodref_base = bodref_base + ctr->nbod;
	cuxUploadConst("dlog", this->dlog);

	lastongpu = true;
}

void host_eventlog::prepare_for_cpu()
{
	if(!lastongpu) { return; }

	// last write was on the GPU, which moved the ref base.
	// before that, we may have had writes on the GPU
	// so now trigger a flush, CPU first, followed by the GPU
	// After the flush is done, download the GPU refbases
	// and store it here

	flush();
}

void host_eventlog::attach_sink(writer *w_)
{
	w = w_;
}

// aux. function to construct the argument-passing structure for writer::process()
void host_eventlog::flush_to_writer()
{
	assert(w);

	output_buffers ob;

	ob.events = events.evt;
	ob.nevt = std::min(ctr->nevt, ecap);
	ob.nevt_dropped = std::max(ctr->nevt - ecap, 0);

	ob.bodies = bodies;
	ob.nbod = std::min(ctr->nbod, bcap);
	ob.nbod_dropped = std::max(ctr->nbod - bcap, 0);

	w->process(ob);

	evtref_base += ctr->nevt;
	bodref_base += ctr->nbod;
	ctr->reset();
}

void host_eventlog::flush()
{
	if(lastongpu)
	{
		// flush CPU buffers
		flush_to_writer();

		// download GPU buffers
		copyFromGPU();
	}
	else
	{
		// CPU was written-to last -- by construction
		// of this class, there's nothing on the GPU side
	}

	flush_to_writer();

	// clear GPU counters
	memcpyToGPU(dlog.ctr, ctr);

	// cpu is now the master
	lastongpu = false;
}

host_eventlog::host_eventlog()
{
	// set everything to NULL
	memset(this, 0, sizeof(this));
}

//
// Initialize the GPU eventlog buffer, and prepare the CPU
// eventlog object
//
void host_eventlog::initialize(int ecap_, int bcap_, int scap_)
{
	if(ctr != NULL) { return; }

	// buffer sizes
	bcap = bcap_;
	ecap = ecap_;
	btresh = (int)(0.9 * bcap);
	etresh = (int)(0.9 * ecap);

	// CPU side
	ctr = new eventlog::counters;
	ctr->reset();
	events.evt = hostAlloc(events.evt, ecap);
	bodies = hostAlloc(bodies, bcap);

	// GPU side
	dlog = *this;
	dlog.ctr = cuxNew<eventlog::counters>(1);
	memcpyToGPU(dlog.ctr, ctr);

	// initialize GPU buffers & eventlog structure
	dlog.events.evt = cuxNew<event>(ecap);
	dlog.bodies = cuxNew<body>(bcap);

	// upload to const
	cuxUploadConst("dlog", this->dlog);
	lastongpu = true;
}

host_eventlog::~host_eventlog()
{
	// Host
	hostFree(events.evt);
	hostFree(bodies);
	delete ctr;
	
	// GPU
	cudaFree(dlog.events.evt);
	cudaFree(dlog.bodies);
	cudaFree(dlog.ctr);
}

//
// gpuPrintf: CPU side
//

bool ievent::isprintf() const
{
	// TODO: WTF?! This shouldn't even compile (it's declared as a const
	// function, and calls the >> extractor ?!?!?!??!?!?!?!)
	int evtid2;
	if(evtid() == EVT_MSGLOST  && *this >> evtid2 && evtid2 == EVT_PRINTF)
	{
		return true;
	}
	if(evtid() != EVT_PRINTF)
	{
		return false;
	}
	return true;
}

std::string ievent::printf() const
{
	// format of data packets:
	// ([hdr][size|...data...][size|...data...])
	std::string res;

	ievent evt(*this->evt);

	int evtid2;
	if(evtid() == EVT_MSGLOST  && evt >> evtid2 && evtid2 == EVT_PRINTF)
	{
		return "Message lost (too big).\n";
	}
	if(evtid() != EVT_PRINTF)
	{
		return "";
	}

	// slurp up the format string
	char fmtbuf[event::SIZEOF], *fmt = fmtbuf;
	evt >> fmt;

	// Now run through it, printing everything we can. We must
	// run to every % character, extract only that, and use printf
	// to format it.
	char buf[1024];
	std::ostringstream out;
	char *p = strchr ( fmt, '%' );
	while ( p != NULL )
	{
		// Print up to the % character
		*p = '\0';
		out << fmt;
		*p = '%';           // Put back the %

		// Now handle the format specifier
		char *format = p++;         // Points to the '%'
		p += strcspn ( p, "%cdiouxXeEfgGaAnps" );
		if ( *p == '\0' )           // If no format specifier, print the whole thing
		{
			fmt = format;
			break;
		}

		ievent::opaque_data m;
		evt >> m;

		char specifier = *p++;
		char c = *p;        // Store for later
		*p = '\0';
		switch ( specifier )
		{
				// These all take integer arguments
			case 'c':
			case 'd':
			case 'i':
			case 'o':
			case 'u':
			case 'x':
			case 'X':
			case 'p':
				sprintf(buf, format, *(int*)m.data);
				out << buf;
				break;

				// These all take double arguments
			case 'e':
			case 'E':
			case 'f':
			case 'g':
			case 'G':
			case 'a':
			case 'A':
				if ( m.size == 4 )  // Float vs. Double thing
					sprintf(buf, format, * ( ( float * ) m.data ));
				else
					sprintf(buf, format, * ( ( double * ) m.data ));
				out << buf;
				break;

				// Strings are handled in a special way
			case 's':
				sprintf(buf, format, ( char * ) m.data);
				out << buf;
				break;

				// % is special
			case '%':
				out << "%%";
				break;

				// Everything else is just printed out as-is
			default:
				out << format;
				break;
		}
		*p = c;                     // Restore what we removed
		fmt = p;                    // Adjust fmt string to be past the specifier
		p = strchr(fmt,'%');    // and get the next specifier
	}

	// Print out the last of the string
	out << fmt;
	return out.str();
}

#if 0
bool next_message(ieventstream &evt, std::string &res)
{
	while(!get_message(evt, res) && evt.next());
	return evt;
}

#endif

//
// Default binary writer
//

class binary_writer : public writer
{
protected:
	std::auto_ptr<std::ostream> eout_strm, bout_strm;
	std::auto_ptr<obstream> eout, bout;
	int ctr;

public:
	binary_writer(const std::string &cfg);
	virtual void process(const output_buffers &ob);
};

extern "C" writer *create_writer_binary(const std::string &cfg)
{
	return new binary_writer(cfg);
}

binary_writer::binary_writer(const std::string &cfg)
{
	ctr = 0;

	std::string event_fn, bodies_fn;
	std::istringstream ss(cfg);
	if(!(ss >> event_fn >> bodies_fn))
		ERROR("Expected 'binary <events_fn> <bodies_fn>' form of configuration for writer.")

	// TODO: check error return codes
	eout_strm.reset(new std::ofstream(event_fn.c_str()));
	bout_strm.reset(new std::ofstream(bodies_fn.c_str()));
	eout.reset(new obstream(*eout_strm));
	bout.reset(new obstream(*bout_strm));
}

// store the accumulated events and bodies into a file, while
// printing out any printf() events to stdout
BLESS_POD(body);
BLESS_POD(event);

void binary_writer::process(const output_buffers &ob)
{
	// download and dump the log
	// write out events
	std::string msg;
	for(int i = 0; i != ob.nevt; i++)
	{
		ievent es(ob.events[i]);
		if(es.isprintf())
		{
			std::string msg = es.printf();
			std::cout << "[Event #" << es.evtref() << " from thread " << es.threadId() << "]: " << msg << "\n";
		}
		else
		{
			// Write it to file
			const event &evt = ob.events[i];
			*eout << evt;

			std::cout << "[Event #" << evt.evtref() << " from thread " << evt.threadId() << "] EVT=" << evt.evtid() << ", " << evt.len() << " bytes long data payload.\n";
		}
	}
	if(ob.nevt_dropped)
	{
		std::cerr << "==== Ran out of event GPU output buffer space (" << ob.nevt_dropped << " event records dropped).\n";
	}

	// write out bodies
	for(int bod = 0; bod != ob.nbod; bod++)
	{
		*bout << ob.bodies[bod];
	}
	if(ob.nbod_dropped)
	{
		std::cerr << "==== Ran out of body GPU output buffer space (" << ob.nbod_dropped << " body records dropped).\n";
	}
}
