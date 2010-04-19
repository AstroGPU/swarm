#include "swarm.h"
#include "swarmio.h"
#include "swarmlog.h"
#include <fstream>
#include <astro/memorymap.h>

// store the accumulated events and bodies into a file, while
// printing out any printf() events to stdout
BLESS_POD(swarm::body);
BLESS_POD(swarm::event);

// aux class to sort indices by energy error (biggest first)
struct energy_sorter
{
	const std::valarray<double> &dE;

	energy_sorter(const std::valarray<double> &dE_) : dE(dE_) {};
	bool operator()(const int  &a, const int  &b) const
	{
		double dEa = fabs(dE[a]);
		double dEb = fabs(dE[b]);
		return dEa > dEb;
	}
};

struct ibfstream
{
	std::ifstream in;
	ibstream bin;

	ibfstream(const std::string &fn)
		: in(fn.c_str()), bin(in)
	{
	}

	template<typename T>
	ibfstream& operator>> (T &v)
	{
		bin >> v;
		return *this;
	}

	void seekg(size_t offs)
	{
		bin.seekg(offs);
	}

	operator bool() const { return bin; }
};

std::ostream &operator<<(std::ostream &out, const swarm::body &b)
{
	char buf[1000];
	sprintf(buf, "%d  % .5e  % .5e % .5e % .5e  % .5e % .5e % .5e\n", b.bod, b.m, b.x, b.y, b.z, b.vx, b.vy, b.vz);
	return out << buf;
}

struct obfstream
{
	const bool binary;
	const bool stdout;
	std::ofstream out;
	obstream bout;

	obfstream(const std::string &fn, bool binary_ = true)
		: stdout(fn == "-"), out(stdout ? "" : fn.c_str()), bout(stdout ? std::cout : out), binary(binary_)
	{
		std::cerr << "Opening " << fn << " for writing.\n";
	}

	template<typename T>
	obfstream& operator<< (const T &v)
	{
		if(binary)
			bout << v;
		else
			(stdout ? std::cout : out) << v;
		return *this;
	}

	operator bool() const { return binary ? (bool)bout : (stdout ? (bool)std::cout : (bool)out); }
};

/*
class log_streamer
{
public:
	virtual void process(body *b, event *evt) = 0;
public:
};
*/

#if 0
/*
	Given the events.bin and bodies.bin file, extract all snapshots (EVT_SNAPSHOT)
	of the system <sys> and store it into the output file <output> in binary or
	text format depending on the value of argument <binary>
*/
int extract_system_snapshots(const std::string &events_fn, const std::string &bodies_fn, int sys, const std::string &output, bool binary, int evtstart = 0)
{
	using namespace swarm;

	ibfstream events(events_fn);
	obfstream out(output, binary);

	// memory-map the bodies
	MemoryMapVector<body> bodies;
	bodies.open(bodies_fn);

	events.seekg(sizeof(event) * evtstart);

	int evtcnt = -1 + evtstart, nsnap = 0;
	event evt;
	while(events >> evt)
	{
		// advance to next snapshot
		evtcnt++;
		if(evt.evtid() != EVT_SNAPSHOT) { continue; }

		// check if this is us
		ievent ie(evt);
		int evsys;
		ie >> evsys;
		if(sys != evsys) { continue; }

		// get system information
		int nbod, snapid, bod_first, bod_last; double T;
		ie >> nbod >> T >> snapid >> bod_first >> bod_last;

		// copy to output
		for(int i = 0, ptr = bod_first; i != nbod; i++)
		{
			body &bod = bodies[ptr];

			// verify consistency
			if(bod.evtref != snapid)
			{
				ERROR("Error: inconsistent events and bodies files around event record " + str(evtcnt) + " and body record number " + str(ptr));
			}

			// copy to output file
			out << bod;

			// move on to the next one
			while(++ptr < bod_last && bodies[ptr].evtref != snapid);
		}

		nsnap++;
	}

	return nsnap;
}

int split_into_per_system(const std::string &events_fn, const std::string &bodies_fn, const std::string &output_fmt, bool binary)
{
	using namespace swarm;

	ibfstream events(events_fn);

	event evt;
	int evtcnt = -1;
	std::set<int> sys_seen;
	while(events >> evt)
	{
		// advance to next snapshot
		evtcnt++;
		if(evt.evtid() != EVT_SNAPSHOT) { continue; }

		// check if we've extracted this one
		ievent ie(evt);
		int sys;
		ie >> sys;
		if(sys_seen.count(sys)) { continue; }

		// extract
		char *fn_tmp;
		asprintf(&fn_tmp, output_fmt.c_str(), sys);
		std::string fn = fn_tmp;
		free(fn_tmp);

		std::cerr << "System " << sys << " -> " << fn << "\n";
		extract_system_snapshots(events_fn, bodies_fn, sys, fn, binary, evtcnt);

		sys_seen.insert(sys);
	}

	return sys_seen.size();
}

#endif

// aux class to sort indices by time (smallest first)
namespace swarm {
struct time_sorter
{

	ensemble &ens;

	time_sorter(ensemble &ens_) : ens(ens_) {};
	bool operator()(const int  &a, const int  &b) const
	{
		return ens.time(a) < ens.time(b);
	}
};
}

// just a simple dump to stdout of one system
void write_output(const swarm::cpu_ensemble &ens, const int sys, std::valarray<double>  &Eold, std::valarray<double> &Enew)
{
	for (int bod = 0; bod != ens.nbod(); bod++)
	{
		float m; double x[3], v[3];
		ens.get_body(sys, bod, m, x[0], x[1], x[2], v[0], v[1], v[2]);

		double dEoverE = (Enew[sys] - Eold[sys]) / Eold[sys];
		printf("%5d %5d  T=%f  m=%f  pos=(% 9.5f % 9.5f % 9.5f)  vel=(% 9.5f % 9.5f % 9.5f)  E=%g  dE/E=%g\n", sys, bod, ens.time(sys), m, x[0], x[1], x[2], v[0], v[1], v[2], Enew[sys], dEoverE);
	}
}

#define SWARMDB 1

int main()
{
	//split_into_per_system("events.bin", "bodies.bin", "system.%04d.txt", false);
	//extract_system_snapshots("events.bin", "bodies.bin", 0, "-", false);
	//return 0;

	using namespace swarm;

	// load the ensemble
	cpu_ensemble ens;

#if SWARMDB
	swarmdb in("log.bin");
	swarmdb::snapshots snaps = in.get_snapshots(swarmdb::ALL);
	snaps.next(ens);
#else
	ens_reader in("output.bin");
	in >> ens;
#endif
	unsigned int nprint = std::min(2, ens.nsys());

	// Calculate energy at beginning of integration
	std::valarray<double> Einit(ens.nsys()), Efinal(ens.nsys());
	calc_total_energy(ens, Einit);

	printf("Snapshot #0 (initial conditions):\n");
	for (unsigned int i = 0;i < nprint;++i)
		write_output(ens, i, Einit, Einit);
	printf("\n");

	// find the final snapshot
#if SWARMDB
	int cnt = 0;
	while(snaps.next(ens)) { cnt++; }
#else
	int cnt = 0;
	while(in >> ens) { cnt++; }
#endif

	// Calculate energy at end of integration
	calc_total_energy(ens, Efinal);

	// store output
	printf("Snapshot #%d (end of simulation)\n", cnt);
	for (unsigned int i = 0;i < nprint;++i)
		write_output(ens, i, Einit, Efinal);
	printf("\n");

	// find systems with worst E conservation
	std::valarray<double> dEoverE = Efinal / Einit - 1.;
	std::vector<int > idx; idx.reserve(ens.nsys());
	for (int i = 0; i != ens.nsys(); i++) idx.push_back(i);
	std::sort(idx.begin(), idx.end(), energy_sorter(dEoverE));
	printf("Systems with worst energy conservation:\n");
	for (unsigned int i = 0;i < nprint;++i)
		write_output(ens, idx[i], Einit, Efinal);

	// find systems with smallest end-of-integration time
	std::sort(idx.begin(), idx.end(), time_sorter(ens));
	printf("\nSystems with smallest time:\n");
	for (unsigned int i = 0;i < nprint;++i)
		write_output(ens, idx[i], Einit, Efinal);

	return 0;
}
