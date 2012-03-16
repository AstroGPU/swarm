//
// C++ Implementation: swarmio
//
// Description: 
//
//
// Author: Mario Juric <mjuric@cfa.harvard.edu>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

/*! \file io.cpp
 *  \brief declares swarm_header
 *
*/

#include "../common.hpp"
#include "io.hpp"


swarm::range_special swarm::ALL;
swarm::range_MIN swarm::MIN;
swarm::range_MAX swarm::MAX;

namespace swarm
{

	/* @file io.cpp
	 *
	 * This file contains routines for opening a swarm log file.
	 * Swarm log file is a binary file with a simple textual header
	 * and a number of fixed size C structs (gpulog::logrecord).
	 *
	 * swarmdb is used to open the log file and query it. API users
	 * should only interact with swarmdb. 
	 *
	 * swarmdb uses indexes for fast retrieval of data. The indexes
	 * are built the first time file is opened and then cached on
	 * disk. There are two indexes:
	 *   1. Time index sorted based on time of records
	 *   2. System index sorted based on system id of records
	 * 
	 * Although sort_binary_output_file function can be used to sort
	 * the entire data file, there is no reason to do so. Since all
	 * the accesses to the file go through the index. 
	 * 
	 *
	 *
	 *
	 *
	 */


	const char* UNSORTED_HEADER_FULL = "unsorted_output // Unsorted output file";
	const char* UNSORTED_HEADER_CHECK = "unsorted_output";
	const char* SORTED_HEADER_FULL = "T_sorted_output // Output file sorted by time";
	const char* SORTED_HEADER_CHECK = "T_sorted_output";
	const char* T_INDEX_CHECK = "T_sorted_index";
	const char* SYS_INDEX_CHECK = "sys_sorted_index";

	void get_Tsys(gpulog::logrecord &lr, double &T, int &sys)
	{
		//std::cerr << "msgid=" << lr.msgid() << "\n";
		if(lr.msgid() < 0)
		{
			// system-defined events that have no (T,sys) heading
			T = -1; sys = -1;
		}
		else
		{
			lr >> T >> sys;
		}
	}

	// Sort the raw outputs
	struct idx_t
	{
		const char *ptr;	// pointer to this packet in memory
		int len;		// length of this packet (including header and data)

		void gethdr(double &T, int &sys) const
		{
			gpulog::logrecord lr(ptr);
			get_Tsys(lr, T, sys);
		}

		bool operator< (const idx_t &a) const
		{
			double Tt, Ta;
			int syst, sysa;
			this->gethdr(Tt, syst);
			a.gethdr(Ta, sysa);

			return 	 Tt < Ta || (
					Tt == Ta && syst < sysa || (
						syst == sysa && this->ptr < a.ptr
					)
				);
		}
	};


 	struct index_entry_time_cmp
 	{
 		bool operator()(const swarmdb::index_entry &a, const swarmdb::index_entry &b) const { return a.T < b.T || (a.T == b.T && a.sys < b.sys || (a.sys == b.sys && a.offs < b.offs) ); }
 	};

 	struct index_entry_sys_cmp
 	{
 		bool operator()(const swarmdb::index_entry &a, const swarmdb::index_entry &b) const { return a.sys < b.sys || (a.sys == b.sys && a.T < b.T) || (a.T == b.T && a.offs < b.offs); }
 	};

#if 0
 	struct index_entry_bod_cmp
 	{
	  bool operator()(const swarmdb::index_entry &a, const swarmdb::index_entry &b) const       { return a.body < b.body || (a.sys == b.sys && a.body ==  b.body && a.T < b.T ) || (a.T == b.T && a.sys == b.sys && a.offs < b.offs); }
 	};
#endif

	bool get_file_info(uint64_t &timestamp, uint64_t &filesize, const std::string &fn)
	{
		struct stat sb;
		if(stat(fn.c_str(), &sb) == -1)
		{
			ERROR(strerror(errno));
		}
#ifdef _WIN32
		timestamp = sb.st_mtime;
#else
		timestamp = (uint64_t(sb.st_mtim.tv_sec) << 32) + uint64_t(sb.st_mtim.tv_nsec);
#endif
		filesize = sb.st_size;
		
		return true;
	}

	template<typename Cmp>
	class index_creator : public index_creator_base
	{
	protected:
		std::string suffix, filetype;
		std::ofstream out;
		std::string filename;
		int nentries;

	public:
		index_creator(const std::string &suffix_, const std::string &filetype_) : suffix(suffix_), filetype(filetype_), nentries(0) {}

		virtual bool start(const std::string &datafile)
		{
			filename = datafile + suffix;
			out.open(filename.c_str());
			assert(out);

			// get the timestamp and file size of the data file
			uint64_t timestamp, filesize;
			get_file_info(timestamp, filesize, datafile);

			// write header
			swarm::swarm_index_header fh(filetype, timestamp, filesize);
			out.write((char*)&fh, sizeof(fh));
			return true;
		}
		virtual bool add_entry(uint64_t offs, gpulog::logrecord lr)
		{
			swarmdb::index_entry ie;
			ie.offs = offs;
			get_Tsys(lr, ie.T, ie.sys);
			out.write((const char *)&ie, sizeof(ie));

			nentries++;
			return true;
		}
		virtual bool finish()
		{
			out.close();

			// sort by time
			mmapped_swarm_index_file mm(filename, filetype, MemoryMap::rw);
			swarmdb::index_entry *begin = (swarmdb::index_entry *)mm.data(), *end = begin + mm.size()/sizeof(swarmdb::index_entry);
			assert((end - begin) == nentries);

			std::sort(begin, end, Cmp());
		#if 1
			Cmp cmp;
			for(int i=1; i < nentries; i++)
			{
				bool ok = cmp(begin[i-1], begin[i]);		// they're less
				if(!ok) { ok = !cmp(begin[i], begin[i-1]); }	// they're equal
				assert(ok);
			}
		#endif
			return true;
		}
	};

	void swarmdb::index_binary_log_file(std::vector<boost::shared_ptr<index_creator_base> > &ic, const std::string &datafile)
	{
		// open datafile
		gpulog::ilogstream ils(mmdata.data(), mmdata.size());
		gpulog::logrecord lr;

		// postprocess (this is where the creator may sort or store the index)
		for(int i=0; i != ic.size(); i++)
		{
			ic[i]->start(datafile);
		}

		// stream through the data file
		while(lr = ils.next())
		{
			for(int i=0; i != ic.size(); i++)
			{
				ic[i]->add_entry(lr.ptr - mmdata.data(), lr);
			}
		}

		// postprocess (this is where the creator may sort or store the index)
		for(int i=0; i != ic.size(); i++)
		{
			ic[i]->finish();
		}
	}

	/*
		Implementation note: This implementation will probably barf
		when log sizes reach a few GB. A better implementation, using merge
		sort could/should be written.
	*/
	bool sort_binary_log_file(const std::string &outfn, const std::string &infn)
	{
		mmapped_swarm_file mm(infn, UNSORTED_HEADER_CHECK);
		gpulog::ilogstream ils(mm.data(), mm.size());
		std::vector<idx_t> idx;
		gpulog::logrecord lr;

		// load record information
		uint64_t datalen = 0;
		for(int i=0; lr = ils.next(); i++)
		{
			idx_t ii;
			ii.ptr = lr.ptr;
			ii.len = lr.len();
			idx.push_back(ii);

			datalen += ii.len;
		}
		assert(datalen == mm.size());

		// sort
		std::sort(idx.begin(), idx.end());

		// output file header
		std::ofstream out(outfn.c_str());
		swarm_header fh(SORTED_HEADER_FULL, 0, datalen);
		out.write((char*)&fh, sizeof(fh));

		// write out the data
		for(int i = 0; i != idx.size(); i++)
		{
			out.write(idx[i].ptr, idx[i].len);
		}
		size_t tp = out.tellp();
		assert(tp == sizeof(fh) + datalen);
		out.close();

		return true;
	}

	struct sysinfo
	{
		int sys;
		long flags;
		const body *bodies;

		bool operator <(const sysinfo &a) const
		{
			return sys < a.sys;
		}
	};

	bool swarmdb::snapshots::next(cpu_ensemble &ens)
	{	
		// find and load the next snapshot
		int nbod = -1, sysmax = -1; double Tsnapend;
		gpulog::logrecord lr;
		std::set<sysinfo> systems;
		while(lr = r.next())
		{
			// find next EVT_SNAPSHOT
			if(lr.msgid() != log::EVT_SNAPSHOT) { continue; }

			// get time & system ID
			double T; int nbod_tmp; sysinfo si;
			lr >> T >> si.sys >> si.flags >> nbod_tmp;
			if(nbod == -1)
			{
				Tsnapend = T*(1+Trelerr) + Tabserr;
				nbod = nbod_tmp;
			}
			assert(nbod == nbod_tmp);	// all systems must have the same number of planets

			// reached the end?
			if(T > Tsnapend)
			{
				r.unget();
				break;
			}
			
			// check if we already have a record of this system
			if(systems.count(si))
			{
				continue;
			}

			// load the pointer to bodies
			lr >> si.bodies;
			systems.insert(si);

			sysmax = std::max(si.sys, sysmax);
		}
		
		// return immediately if we've reached the end
		if(!systems.size()) { return false; }

		// pack everything to cpu_ensemble structure
		ens = cpu_ensemble::create(nbod, sysmax+1);
		// initially, mark everything as inactive
		for(int sys = 0; sys != ens.nsys(); sys++)
		{
			ens.set_inactive(sys);
		}
		
		// populate the systems with data and mark them active
		for(std::set<sysinfo>::iterator i = systems.begin(); i != systems.end(); i++)
		{
			const sysinfo &si = *i;
			assert(si.sys >= 0 && si.sys < ens.nsys());

			// per-system data
			ens.flags(si.sys) = si.flags;

			// per-body data
			for(int bod = 0; bod != nbod; bod++)
			{
				const body &b = si.bodies[bod];
				ens.set_body(si.sys, bod,  b.mass, b.x, b.y, b.z, b.vx, b.vy, b.vz);
				ens.time(si.sys) = Tsnapend;
			}
		}

		return true;
	}

	swarmdb::snapshots::snapshots(const swarmdb &db_, time_range_t T, double Tabserr_, double Trelerr_)
		: db(db_), Tabserr(Tabserr_), Trelerr(Trelerr_), r(db.query(ALL, T))
	{
	}

	swarmdb::result::result(const swarmdb &db_, const sys_range_t &sys_, const time_range_t &T_)
	  : db(db_), sys(sys_), T(T_)
	{
		using namespace boost;

		if(sys.first == sys.last || !T)
		{
			// find the sys index range
			swarmdb::index_entry dummy;
			dummy.sys = sys.first; begin = std::lower_bound(db.idx_sys.begin, db.idx_sys.end, dummy, bind( &index_entry::sys, _1 ) < bind( &index_entry::sys, _2 ));
			dummy.sys = sys.last;  end   = std::upper_bound(db.idx_sys.begin, db.idx_sys.end, dummy, bind( &index_entry::sys, _1 ) < bind( &index_entry::sys, _2 ));
		}
		else if(T)
		{
			// find the first T in the range
			swarmdb::index_entry dummy;
			dummy.T = T.first; begin = std::lower_bound(db.idx_time.begin, db.idx_time.end, dummy, bind( &index_entry::T, _1 ) < bind( &index_entry::T, _2 ));
			dummy.T = T.last;  end   = std::upper_bound(db.idx_time.begin, db.idx_time.end, dummy, bind( &index_entry::T, _1 ) < bind( &index_entry::T, _2 ));
		}
		else
		{
			// stream through everything, time first
			begin = db.idx_time.begin;
			end   = db.idx_time.end;
		}

		at = begin;
		atprev = at;
	}


#if 0
  	swarmdb::result::result(const swarmdb &db_, const sys_range_t &sys_, const body_range_t &body_, const time_range_t &T_)	  : db(db_), sys(sys_), body(body_), T(T_)
	{
		using namespace boost;

		if(sys.first == sys.last || !T)
		{
			// find the sys index range
			swarmdb::index_entry dummy;
			dummy.sys = sys.first; begin = std::lower_bound(db.idx_sys.begin, db.idx_sys.end, dummy, bind( &index_entry::sys, _1 ) < bind( &index_entry::sys, _2 ));
			dummy.sys = sys.last;  end   = std::upper_bound(db.idx_sys.begin, db.idx_sys.end, dummy, bind( &index_entry::sys, _1 ) < bind( &index_entry::sys, _2 ));
		}
		else if(T)
		{
			// find the first T in the range
			swarmdb::index_entry dummy;
			dummy.T = T.first; begin = std::lower_bound(db.idx_time.begin, db.idx_time.end, dummy, bind( &index_entry::T, _1 ) < bind( &index_entry::T, _2 ));
			dummy.T = T.last;  end   = std::upper_bound(db.idx_time.begin, db.idx_time.end, dummy, bind( &index_entry::T, _1 ) < bind( &index_entry::T, _2 ));
		}
		else
		{
			// stream through everything, time first
			begin = db.idx_time.begin;
			end   = db.idx_time.end;
		}

		at = begin;
		atprev = at;
	}
#endif

	gpulog::logrecord swarmdb::result::next()
	{
		if(at < end)
		{
			atprev = at;
		}

		while(at < end)
		{
			if(!T.in(at->T))       { at++; continue; }
			if(!sys.in(at->sys))   { at++; continue; }
			//			if(!body.in(at->body)) { at++; continue; }

			return gpulog::logrecord(db.mmdata.data() + (at++)->offs);
		}

		static gpulog::header hend(-1, 0);
		static gpulog::logrecord eof((char *)&hend);
		return eof;
	}

	void swarmdb::result::unget()
	{
		at = atprev;
	}

	swarmdb::swarmdb(const std::string &datafile)
	{
		open(datafile);
	}

	void swarmdb::open(const std::string &datafile)
	{
		if(access(datafile.c_str(), 04) != 0)
		{
			ERROR("Cannot open output file '" + datafile + "'");
		}
		this->datafile = datafile;

		// open the datafile
		mmdata.open(datafile, UNSORTED_HEADER_CHECK);

		// open the indexes
		open_indexes();
	}

	void swarmdb::open_indexes(bool force_recreate)
	{
		std::vector<boost::shared_ptr<index_creator_base> > ic;
		std::string fn;

		// auto-create indices if needed
		bool tidx_open = !force_recreate && open_index(idx_time, datafile, ".time.idx", "T_sorted_index");
		if(!tidx_open)
		{
			boost::shared_ptr<index_creator_base> ii(new index_creator<index_entry_time_cmp>(".time.idx", "T_sorted_index"));
			ic.push_back( ii );
		}
		bool sysidx_open = !force_recreate && open_index(idx_sys,  datafile, ".sys.idx", "sys_sorted_index");
		if(!sysidx_open)
		{
			boost::shared_ptr<index_creator_base> ii(new index_creator<index_entry_sys_cmp>(".sys.idx", "sys_sorted_index"));
			ic.push_back( ii );
		}

		// create indices if needed
		if(!ic.empty())
		{
			index_binary_log_file(ic, datafile);
		}

		// open index maps
		if(!tidx_open && !open_index(idx_time, datafile, ".time.idx", T_INDEX_CHECK ))
		{
 			ERROR("Cannot open index file '" + datafile + ".time.idx'");
 		}
		if(!sysidx_open && !open_index(idx_sys,  datafile, ".sys.idx", SYS_INDEX_CHECK))
		{
 			ERROR("Cannot open index file '" + datafile + ".sys.idx'");
 		}
	}

	bool swarmdb::open_index(index_handle &h, const std::string &datafile, const std::string &suffix, const std::string &filetype)
	{
		std::string filename = datafile + suffix;

		// check for existence
		if(access(filename.c_str(), 04) != 0) { return false; }

		

		// check for modification timestamp and size
		h.mm.open(filename.c_str(), filetype);
		uint64_t timestamp, filesize;
		get_file_info(timestamp, filesize, datafile);

		//TODO: this quick fix is only for fake memory mapping
		timestamp = h.mm.hdr().timestamp;

		if(h.mm.hdr().datafile_size != filesize || h.mm.hdr().timestamp != timestamp)
		{
			std::cerr << "Index " << filename << " not up to date. Will regenerate.\n";
			return false;
		}
	
		h.begin = (swarmdb::index_entry *)h.mm.data();
		h.end   = h.begin + h.mm.size()/sizeof(swarmdb::index_entry);
		
		return true;
	}
}
