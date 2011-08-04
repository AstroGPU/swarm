/*! \file io.hpp
 *  \brief declares swarmdb and range
 *
*/

#ifndef swarmio_h__
#define swarmio_h__

#include "swarm.h"
#include "log.hpp"
#include "peyton/binarystream.hpp"
#include "peyton/memorymap.hpp"

#include <fstream>
#include <sstream>
#include <limits>
#include "fileformat.hpp"
namespace swarm {

	extern struct range_special { } ALL;
	extern struct range_MAX
	{
		template<typename T> operator T() const
		{
			return std::numeric_limits<T>::max();
		};
	} MAX;
	extern struct range_MIN
	{
		template<typename T> operator T() const
		{
			return std::numeric_limits<T>::is_integer ? std::numeric_limits<T>::min() : -std::numeric_limits<T>::max();
		};
	} MIN;

	template<typename T>
	struct range
	{
		T first, last;

		range(const T &a) : first(a), last(a) {}
		range(const T &a, const T &b) : first(a), last(b) {}
		range(const range_special &r = ALL) : first(MIN), last(MAX) {}

		bool in(const T& v) { return first <= v && v <= last; }
		operator bool() const { return first <= last; }
	};

	typedef range<int> sys_range_t;
	typedef range<double> time_range_t;


	typedef mmapped_file_with_header<swarm_header> mmapped_swarm_file;
	typedef mmapped_file_with_header<swarm_index_header> mmapped_swarm_index_file;

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
			mmapped_swarm_index_file mm;
			const index_entry *begin, *end;
		};

		mmapped_swarm_file mmdata;

		index_handle idx_time, idx_sys;
		std::string datafile;

		void open(const std::string &datafile);
		void open_indexes(bool force_recreate = false);
		bool open_index(index_handle &h, const std::string &datafile, const std::string &suffix, const std::string &filetype);

	public:
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

			bool next(cpu_ensemble &ens/*, bool keep_existing = true*/);
			snapshots(const swarmdb &db, time_range_t T, double Tabserr = 0, double Trelerr = 0);
		};

		snapshots get_snapshots(time_range_t T, double Tabserr = 0, double Trelerr = 0)
		{
			return snapshots(*this, T, Tabserr, Trelerr);
		}
	};

	bool sort_binary_log_file(const std::string &outfn, const std::string &infn);

} // end namespace swarm

#endif
