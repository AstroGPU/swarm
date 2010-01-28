#include "swarmlog.h"
#include <iostream>
#include <astro/binarystream.h>
#include <fstream>
#include <sstream>
#include <memory>

void debug_hook()
{
	std::cerr << "";
}

cpu_eventlog clog;
static bool cpu_eventlog_initialized = 0;

const ieventstream::body* ieventstream::get_bodies(int &nbod, int &ndropped) const
{
	nbod = std::min(clog.ctr->nbodX, clog.bcapX);
	ndropped = std::max(clog.ctr->nbodX - clog.bcapX, 0);
	return clog.bodies;
}

void initialize_eventlog(int ecap, int bcap, int scap)
{
	clog.initialize(ecap, bcap, scap);
}

ieventstream::ieventstream(bool syncNew)
{
	assert(cpu_eventlog_initialized);

	if(syncNew) { clog.sync(); }
	init(clog.events, clog.ctr->nevtX, clog.ecapX);
}

// advance to next message
int ieventstream::next()
{
	// already at end?
	if(atmsg == nmsg) { return eventlog_base::EVT_EOF; }

	// move and exit if we're at the end
	atmsg++;
	if(atmsg != 0)
	{
		data += eventlog_base::MAX_MSG_LEN;
	}

	// learn the basics about the current message
	if(atmsg != nmsg)
	{
		hdr = *(eventlog_base::evt_hdr*)(data);
		at = sizeof(eventlog_base::evt_hdr);
		m_eom = false;
	}
	else
	{
		hdr.evtid = eventlog_base::EVT_EOF;
	}

	// return current message event ID
	return evtid();
}

void ieventstream::init(const char *data_, int nmsg_, int ecap_)
{
	data = data_;
	nmsg = std::min(nmsg_, ecap_);
	m_nlost = std::max(nmsg_ - ecap_, 0);

	// position us _ahead_ of the first message
	// (expecting a call to next() will move us to first message)
	atmsg = -1;
	m_eom = true;
}

bool ieventstream::check_end()
{
	if(atmsg == -1) // first call, auto-advance
	{
		next();
	}
	if(hdr.len == at) { m_eom = true; }
	return *this;
}

ieventstream &ieventstream::operator >>(raw_msg &v)
{
	if(!check_end()) { return *this; }

	v.size = *(int*)(data + at);
	at += sizeof(int);
	v.data = data+at;
	at += v.size;
	return *this;
}

ieventstream &ieventstream::operator >>(std::string &s)
{
	if(!check_end()) { return *this; }

	int size = *(int*)(data + at);
	at += sizeof(int);
	s = (data+at);
	int strlen = s.size();
	if(strlen != size-1-sizeof(int))
		ERROR("Programmer error: data size != string size. Contact the authors.");
	at += size;
	return *this;
}

ieventstream &ieventstream::operator >>(char *s)
{
	if(!check_end()) { return *this; }

	int size = *(int*)(data + at);
	at += sizeof(int);
	strcpy(s, data+at);
	int strl = strlen(s);
	if(strl != size-1)
		ERROR("Programmer error: data size != string size. Contact the authors.");
	at += size;
	return *this;
}

ieventstream &ieventstream::operator >>(eventlog_base::event &evt)
{
	if(!check_end()) { return *this; }
	assert(sizeof(eventlog_base::event) == eventlog_base::MAX_MSG_LEN);
	memcpy(&evt, data, sizeof(eventlog_base::event));
	at = hdr.len;

	return *this;
}

// download and clear the contents of GPU log
void cpu_eventlog::sync()
{
	// download GPU log data
	if(!cpu_eventlog_initialized)
		ERROR("Programmer error: cpu_eventlog must be initialized before use. Contact the authors.");

	// download GPU counters
	memcpyToHost(ctr, glog.ctr);

	// download the data
	events = hostAlloc(events, ecapX * MAX_MSG_LEN);
	memcpyToHost(events, glog.events, ecapX * MAX_MSG_LEN);
	bodies = hostAlloc(bodies, bcapX);
	memcpyToHost(bodies, glog.bodies, bcapX);

	// clear GPU counters
	counters tmp = *ctr;
	tmp.reset();
	memcpyToGPU(glog.ctr, &tmp);
}

void cpu_eventlog::initialize(int ecap_, int bcap_, int scap_)
{
	if(cpu_eventlog_initialized) { return; }
	cpu_eventlog_initialized = true;

	// buffer sizes
	bcapX = bcap_;
	ecapX = ecap_;
	btresh = (int)(0.9 * bcapX);
	etresh = (int)(0.9 * ecapX);

	// copy self to GPU version
	glog = *this;

	// now update the GPU version with its own pointers
	// initialize counters/pointers
	ctr = new eventlog_base::counters;
	ctr->reset();
	glog.ctr = cuxNew<eventlog_base::counters>(1);
	memcpyToGPU(glog.ctr, ctr);

	// initialize buffers & eventlog structure
	glog.events = cuxNew<char>(ecapX * MAX_MSG_LEN);
	glog.bodies = cuxNew<body>(bcapX);
	//glog.systems = NULL;

	// upload to const
	cuxUploadConst("glog", this->glog);
}


bool get_as_message(ieventstream &msg, std::string &res);

bool next_message(ieventstream &msg, std::string &res)
{
	while(!get_as_message(msg, res) && msg);
	return msg;
}

bool get_as_message(ieventstream &msg, std::string &res)
{
	// format of data packets:
	// ([hdr][size|...data...][size|...data...])

	int evtid2;
	if(msg.evtid() == cpu_eventlog::EVT_MSGLOST  && msg >> evtid2 && evtid2 == cpu_eventlog::EVT_PRINTF)
	{
		msg.next();
		res = "Message lost (too big).\n";
		return true;
	}
	if(msg.evtid() != cpu_eventlog::EVT_PRINTF)
	{
		msg.next();
		return false;
	}

	// slurp up the format string
	char fmtbuf[cpu_eventlog::MAX_MSG_LEN], *fmt = fmtbuf;
	msg >> fmt;

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

		ieventstream::raw_msg m;
		msg >> m;

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
	res = out.str();
	return true;
}

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
	virtual void process(ieventstream &es);
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
BLESS_POD(eventlog_base::body);
void binary_writer::process(ieventstream &es)
{
	// download and dump the log
	// write out events
	while(es.next())
	{
		std::string msg;
		if(get_as_message(es, msg))
		{
			std::cout << "[Message (" << ctr++ << ")]: " << msg << "\n";
		}
		else
		{
			// Write it to file
			std::cout << "[Event " << es.evtid() << "]\n";
			eventlog_base::event evt;
			es >> evt;
			*eout << evt;
		}
	}
	if(es.nevents_dropped())
	{
		std::cerr << "==== Ran out of event GPU output buffer space (" << es.nevents_dropped() << " event records dropped).\n";
	}

	// write out bodies
	int nbod, ndropped;
	const eventlog_base::body *bodies = es.get_bodies(nbod, ndropped);
	for(int bod = 0; bod != nbod; bod++)
	{
		*bout << bodies[bod];
	}
	if(ndropped)
	{
		std::cerr << "==== Ran out of body GPU output buffer space (" << ndropped << " body records dropped).\n";
	}
}

void test()
{
	cpu_eventlog elog;
}
