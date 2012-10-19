/*************************************************************************
 * Copyright (C) 2011 by Eric Ford and the Swarm-NG Development Team     *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 3 of the License.        *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program; if not, write to the                         *
 * Free Software Foundation, Inc.,                                       *
 * 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ************************************************************************/

/*! \file io.hpp
 *  \brief Defines swarmdb and range.
 *
*/

#ifndef swarmio_h__
#define swarmio_h__

#include "../peyton/binarystream.hpp"
#include "../peyton/memorymap.hpp"

#include "fileformat.hpp"
#include "log.hpp"

namespace swarm {

	extern const char* UNSORTED_HEADER_FULL;
	extern const char* UNSORTED_HEADER_CHECK;
	extern const char* SORTED_HEADER_FULL;
	extern const char* SORTED_HEADER_CHECK;

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
	typedef range<int> body_range_t;
	typedef range<double> time_range_t;

	typedef mmapped_file_with_header<swarm_header> mmapped_swarm_file;
	typedef mmapped_file_with_header<swarm_index_header> mmapped_swarm_index_file;

	struct index_creator_base
	{
		virtual bool start(const std::string &datafile) = 0;
		virtual bool add_entry(uint64_t offs, gpulog::logrecord lr) = 0;
		virtual bool finish() = 0;
		virtual ~index_creator_base() {};
	};

	class swarmdb
	{
	public:
		struct index_entry
		{
			uint64_t offs;	// data offset for the record
	
			double T;	// time
			int sys;	// system at the record
			int body;	// bod at/in the record
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
		  //			body_range_t  body;
			time_range_t T;
			
			const index_entry *begin, *end, *at, *atprev;

		  result(const swarmdb &db_, const sys_range_t &sys, const time_range_t &T);
		  //		  result(const swarmdb &db_, const sys_range_t &sys, const body_range_t &body, const time_range_t &T);

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

	  /*
	  result query(sys_range_t sys, body_range_t body, time_range_t T) const
		{
		  return result(*this, sys, body, T);
		}
	  */

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
	private:
		void index_binary_log_file(std::vector<boost::shared_ptr<index_creator_base> > &ic, const std::string &datafile);
	};

	bool sort_binary_log_file(const std::string &outfn, const std::string &infn);

} // end namespace swarm

#endif
