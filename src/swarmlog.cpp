#include "swarmlog.h"
#include <iostream>

void debug_hook()
{
	std::cerr << "";
}

// CPU end of the event log -- receives the events recorded
// on the GPU
struct cpueventlog : public eventlog
{
protected:
	eventlog glog;

public:
	void initialize(int ecap = 16*1024, int bcap = 1024*1024, int scap = 1024*1024);

	void sync();

	friend ieventstream get_gpu_events();
};

static cpueventlog clog;
static bool cpueventlog_initialized = 0;

void initialize_eventlog(int ecap, int bcap, int scap)
{
	clog.initialize(ecap, bcap, scap);
}

ieventstream::ieventstream(bool syncNew)
{
	assert(cpueventlog_initialized);

	if(syncNew) { clog.sync(); }
	init(clog.events, clog.ctr->nevtX, clog.ecapX);
}

// advance to next message
int ieventstream::next()
{
	// already at end?
	if(atmsg == nmsg) { return eventlog::EVT_EOF; }

	// move and exit if we're at the end
	atmsg++;
	if(atmsg != 0)
	{
		data += eventlog::MAX_MSG_LEN;
	}

	// learn the basics about the current message
	if(atmsg != nmsg)
	{
		hdr = *(eventlog::evt_hdr*)(data);
		at = sizeof(eventlog::evt_hdr);
		m_eom = false;
	}
	else
	{
		hdr.evtid = eventlog::EVT_EOF;
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

// download and clear the contents of GPU log
void cpueventlog::sync()
{
	// download GPU log data
	if(!cpueventlog_initialized)
		ERROR("Programmer error: cpueventlog must be initialized before use. Contact the authors.");

	// download GPU counters
	memcpyToHost(ctr, glog.ctr);

	// download the data
	events = hostAlloc(events, ecapX * MAX_MSG_LEN);
	memcpyToHost(events, glog.events, ecapX * MAX_MSG_LEN);

	// reset GPU counters
	counters tmp = *ctr;
	tmp.reset();
	memcpyToGPU(glog.ctr, &tmp);
}

void cpueventlog::initialize(int ecap_, int bcap_, int scap_)
{
	if(cpueventlog_initialized) { return; }
	cpueventlog_initialized = true;

	// buffer sizes
	bcap = scap = 0; ecapX = ecap_;
	btresh = stresh = 0; etresh = (int)(0.9 * ecapX);

	// copy self to GPU version
	glog = *this;

	// now update the GPU version with its own pointers
	// initialize counters/pointers
	ctr = new eventlog::counters;
	ctr->reset();
	glog.ctr = cuxNew<eventlog::counters>(1);
	memcpyToGPU(glog.ctr, ctr);

	// initialize buffers & eventlog structure
	glog.events = cuxNew<char>(ecapX * MAX_MSG_LEN);
	glog.bodies = NULL;
	glog.systems = NULL;

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
	if(msg.evtid() == cpueventlog::EVT_MSGLOST  && msg >> evtid2 && evtid2 == cpueventlog::EVT_PRINTF)
	{
		msg.next();
		res = "Message lost (too big).\n";
		return true;
	}
	if(msg.evtid() != cpueventlog::EVT_PRINTF)
	{
		msg.next();
		return false;
	}

	// slurp up the format string
	char fmtbuf[cpueventlog::MAX_MSG_LEN], *fmt = fmtbuf;
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

void dump_gpu_events()
{
	// download and dump the log
	// NOTE: this should be done less frequently
	ieventstream es;
	static int ctr = 0;
	while(es.next())
	{
		std::string msg;
		if(get_as_message(es, msg))
		{
			std::cout << "[Message (" << ctr++ << ")]: " << msg << "\n";
		}
		else
		{
			std::cout << "[Event " << es.evtid() << "]\n";
		}
	}
	if(es.nlost())
	{
		std::cout << "==== Ran out of GPU event buffer space. " << es.nlost() << " events lost at this point.\n";
	}
	else
	{
		std::cout << "====\n";
	}
}
