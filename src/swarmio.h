#ifndef swarmio_h__
#define swarmio_h__

#include "swarm.h"
#include "swarmlog.h"
#include <astro/binarystream.h>
#include <astro/memorymap.h>
#include <fstream>
#include <sstream>
#include <limits>

// Convert a variable of arbitrary type to a string.
// NOTE: heavy (unoptimized) function, use sparingly
template<typename T>
std::string str(const T& var)
{
	std::ostringstream ss;
	ss << var;
	return ss.str();
}

namespace swarm {

//
// I/O and snapshotting support
//
class ens_writer
{
protected:
	std::string fn;
	std::ofstream out;
	obstream bout;
public:
	ens_writer(const std::string &fn);
	ens_writer &operator <<(const cpu_ensemble &ens);
	operator bool() const { return bout; }
};

class ens_reader
{
protected:
	std::string fn;
	std::ifstream in;
	ibstream bin;
public:
	ens_reader(const std::string &fn);
	ens_reader &operator >>(cpu_ensemble &ens);
	operator bool() const { return bin; }
};

	struct file_header;
	struct mmapped_swarm_file : public MemoryMap
	{
	protected:
		file_header *fh;
		char *m_data;
		int m_size;

	public:
		mmapped_swarm_file(const std::string &filename = "", const std::string &type = "", int mode = ro, bool validate = true);
		void open(const std::string &filename, const std::string &type, int mode = ro, bool validate = true);
		
		// override MemoryMap::size() to return the length of the data without the header
		size_t size() const { return m_size; }

		// accessors
		char *data() const { return m_data; }
		file_header &hdr() const { return *fh; }
	};

	class swarmdb
	{
	public:
		struct index_entry
		{
			uint64_t offs;	// data offset for the record
	
			double T;	// time
			int sys;	// system at the record
		};

	protected:
		struct index_handle
		{
			mmapped_swarm_file mm;
			const index_entry *begin, *end;
		};

		mmapped_swarm_file mmdata;

		index_handle idx_time, idx_sys;
		std::string datafile;

		void open(const std::string &datafile);
		void open_indexes(bool recreate = true);
		void open_index(index_handle &h, const std::string &idxfile, const std::string &filetype);

	public:
		static struct range_special { } ALL;
		static struct range_MAX
		{
			template<typename T> operator T() const
			{
				return std::numeric_limits<T>::max();
			};
		} MAX;
		static struct range_MIN
		{
			template<typename T> operator T() const
			{
				return std::numeric_limits<T>::is_integer ? std::numeric_limits<T>::min() : -std::numeric_limits<T>::max();
			};
		} MIN;

		template<typename T>
		struct range
		{
			T begin, end;

			range(const T &a) : begin(a), end(a + 1) {}
			range(const T &a, const T &b) : begin(a), end(b) {}
			range(const range_special &r) : begin(MIN), end(MAX) {}

			bool in(const T& v) { return begin <= v && v < end; }
			operator bool() const { return begin < end; }
			ptrdiff_t width() const { return end - begin; }
		};

		typedef range<int> sys_range_t;
		typedef range<double> time_range_t;

		struct result
		{
			const swarmdb &db;

			sys_range_t  sys;
			time_range_t T;
			
			const index_entry *begin, *end, *at, *atprev;

			result(const swarmdb &db_, const sys_range_t &sys, const time_range_t &T);

			gpulog::logrecord next();
			void unget();
		};

	public:
		swarmdb(const std::string &datafile);

		// return a stream of events with msgid, and system sys, at time T
		result query(sys_range_t sys, time_range_t T) const
		{
			return result(*this, sys, T);
		}

	public:
		struct snapshots
		{
			const swarmdb &db;
			result r;

			double Tabserr, Trelerr;

			bool next(cpu_ensemble &ens, bool keep_existing = true);
			snapshots(const swarmdb &db, time_range_t T, double Tabserr = 0, double Trelerr = 0);
		};

		snapshots get_snapshots(time_range_t T, double Tabserr = 0, double Trelerr = 0)
		{
			return snapshots(*this, T, Tabserr, Trelerr);
		}
	};


} // end namespace swarm

#endif
