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

#include "swarmlog.h"
#include "swarmio.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <limits>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace swarm
{
	// Note: the header _MUST_ be padded to 16-byte boundary
	struct ALIGN(16) file_header
	{
		char magic[6];		// Magic string to quickly verify this is a swarm file (== 'SWARM\0')
		char version[2];	// File format version
		char m_type[76];	// user-defined file content ID/description
		uint32_t flags;		// user-defined flags
		uint64_t datalen;	// length of data in the file (0xFFFFFFFFFFFFFFFF for unknown)

		static const uint64_t npos = 0xFFFFFFFFFFFFFFFFLL;

		file_header(const std::string &type, int flags_ = 0, uint64_t datalen_ = npos)
		{
			strcpy(magic, "SWARM");
			strcpy(version, "0");

			strncpy(this->m_type, type.c_str(), sizeof(this->m_type));
			this->m_type[sizeof(this->m_type)-1] = '\0';

			flags = flags_;
			datalen = datalen_;
		}

		std::string type() const
		{
			char *c = strstr(m_type, "//");
			if(!c) { return trim(m_type); }
			return trim(std::string(m_type, c - m_type));
		}

		bool is_compatible(const file_header &a)
		{
			bool ffver_ok = memcmp(magic, a.magic, 8) == 0;
			std::string t = type();
			bool type_ok = t.empty() || (t == a.type());
			return ffver_ok && type_ok;
		}
	};

	mmapped_swarm_file::mmapped_swarm_file(const std::string &filename, const std::string &type, int mode, bool validate)
		: fh(NULL), m_data(NULL)
	{
		if(!filename.empty())
		{
			open(filename, type, mode, validate);
		}
	}

	void mmapped_swarm_file::open(const std::string &filename, const std::string &type, int mode, bool validate)
	{
		MemoryMap::open(filename.c_str(), 0, 0, mode, shared);

		// read/check header
		fh = (file_header *)map;
		file_header ref_ver(type);
		if(validate && !ref_ver.is_compatible(*fh))
		{
			std::cerr << "Expecting: " << ref_ver.type() << "\n";
			std::cerr << "Got      : " << fh->type() << "\n";
			ERROR("Input file corrupted or incompatible with this version of swarm");
		}

		// get the pointer to data
		m_data = (char *)&fh[1];
		
		// data size
		m_size = length - sizeof(file_header);
	}

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
		uint64_t offs;		// offset to this packet in memory
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

			return 	 Tt < Ta ||
				(Tt == Ta && syst < sysa);
		}
	};

	struct index_creator_base
	{
		virtual bool start(const std::string &datafile) = 0;
		virtual bool add_entry(uint64_t offs, gpulog::logrecord &lr) = 0;
		virtual bool finish() = 0;
		virtual ~index_creator_base() {};
	};

 	struct index_entry_time_cmp
 	{
 		bool operator()(const swarmdb::index_entry &a, const swarmdb::index_entry &b) const { return a.T < b.T || (a.T == b.T && a.sys < b.sys || (a.sys == b.sys && a.offs < b.offs) ); }
 	};

 	struct index_entry_sys_cmp
 	{
 		bool operator()(const swarmdb::index_entry &a, const swarmdb::index_entry &b) const { return a.sys < b.sys || (a.sys == b.sys && a.T < b.T) || (a.T == b.T && a.offs < b.offs); }
 	};

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

			// write header
			swarm::file_header fh(filetype);
			out.write((char*)&fh, sizeof(fh));
		}
		virtual bool add_entry(uint64_t offs, gpulog::logrecord &lr)
		{
			swarmdb::index_entry ie;
			ie.offs = offs;
			get_Tsys(lr, ie.T, ie.sys);
			out.write((const char *)&ie, sizeof(ie));
			nentries++;
		}
		virtual bool finish()
		{
			out.close();

			// sort by time
			mmapped_swarm_file mm(filename, filetype, MemoryMap::rw);
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
		}
	};

	void index_binary_log_file(std::vector<boost::shared_ptr<index_creator_base> > &ic, const std::string &datafile)
	{
		// open datafile
		mmapped_swarm_file mm(datafile, "T_sorted_output");
		gpulog::ilogstream ils(mm.data(), mm.size());
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
				ic[i]->add_entry(lr.ptr - mm.data(), lr);
			}
		}

		// postprocess (this is where the creator may sort or store the index)
		for(int i=0; i != ic.size(); i++)
		{
			ic[i]->finish();
		}
	}

	swarmdb::range_special swarmdb::ALL;
	swarmdb::range_MIN swarmdb::MIN;
	swarmdb::range_MAX swarmdb::MAX;

	/*
		Implementation note: This implementation will probably barf
		when log sizes reach a few GB. A better implementation, using merge
		sort could/should be written.
	*/
	bool sort_binary_log_file(const std::string &outfn, const std::string &infn)
	{
		mmapped_swarm_file mm(infn, "unsorted_output");
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
		file_header fh("T_sorted_output // Output file sorted by time", 0, datalen);
		out.write((char*)&fh, sizeof(fh));

		// write out the data
		uint64_t at = 0;
		for(int i = 0; i != idx.size(); i++)
		{
			out.write(idx[i].ptr, idx[i].len);
			idx[i].offs = at;
		}
		size_t tp = out.tellp();
		assert(tp == sizeof(fh) + datalen);
		out.close();

		return true;
	}

	struct sysinfo
	{
		int sys, flags;
		const body *bodies;
	};

	bool swarmdb::snapshots::next(cpu_ensemble &ens, bool keep_existing)
	{	
		// find and load the next snapshot
		int nbod = -1, sysmax = -1; double Tsnapend;
		gpulog::logrecord lr;
		std::vector<sysinfo> systems;
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

			// load the pointer to bodies
			lr >> si.bodies;
			systems.push_back(si);

			sysmax = std::max(si.sys, sysmax);
		}
		
		// return immediately if we've reached the end
		if(!systems.size()) { return false; }

		// pack everything to cpu_ensemble structure
		if(keep_existing && ens.nsys() != 0)
		{
			assert(sysmax < ens.nsys());
		}
		else
		{
			ens.reset(sysmax+1, nbod);
		}

		// initially, mark everything as inactive
		for(int sys = 0; sys != ens.nsys(); sys++)
		{
			ens.set_inactive(sys);
		}
		
		// populate the systems with data and mark them active
		for(int i = 0; i != systems.size(); i++)
		{
			sysinfo &si = systems[i];
			assert(si.sys >= 0 && si.sys < ens.nsys());

			// per-system data
			ens.flags(si.sys) = si.flags;

			// per-body data
			for(int bod = 0; bod != nbod; bod++)
			{
				const body &b = si.bodies[bod];
				ens.set_body(si.sys, bod,  b.m, b.x, b.y, b.z, b.vx, b.vy, b.vz);
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
		if(sys.width() == 1 || !T)
		{
			// find the sys index range
			swarmdb::index_entry dummy = { 0, 0, 0 };
			dummy.sys = sys,begin; begin = std::lower_bound(db.idx_sys.begin, db.idx_sys.end, dummy, index_entry_sys_cmp());
			dummy.sys = sys.end;   end   = std::upper_bound(db.idx_sys.begin, db.idx_sys.end, dummy, index_entry_sys_cmp());
		}
		else if(T)
		{
			// find the first T in the range
			swarmdb::index_entry dummy = { 0, 0, 0 };
			dummy.T = T.begin; begin = std::lower_bound(db.idx_time.begin, db.idx_time.end, dummy, index_entry_time_cmp());
			dummy.T = T.end;   end   = std::upper_bound(db.idx_time.begin, db.idx_time.end, dummy, index_entry_time_cmp());
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

	gpulog::logrecord swarmdb::result::next()
	{
		if(at < end)
		{
			atprev = at;
		}

		while(at < end)
		{
			if(T   && !T.in(at->T))     { at++; continue; }
			if(sys && !sys.in(at->sys)) { at++; continue; }

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
		if(access(datafile.c_str(), R_OK) != 0)
		{
			ERROR("Cannot open output file '" + datafile + "'");
		}
		this->datafile = datafile;

		// open the datafile
		mmdata.open(datafile, "T_sorted_output");

		// open the indexes
		open_indexes();
	}

	void swarmdb::open_indexes(bool recreate)
	{
		std::vector<boost::shared_ptr<index_creator_base> > ic;
		std::string fn;

		// auto-create indices if needed
		if(recreate || access((datafile + ".time.idx").c_str(), R_OK) != 0)
		{
			ic.push_back( boost::make_shared< index_creator<index_entry_time_cmp> >(".time.idx", "T_sorted_index") );
		}
		if(recreate || access((datafile + ".sys.idx").c_str(), R_OK) != 0)
		{
			ic.push_back( boost::make_shared< index_creator<index_entry_sys_cmp> >(".sys.idx", "sys_sorted_index") );
		}
		if(!ic.empty())
		{
			index_binary_log_file(ic, datafile);
		}

		// open index maps
		open_index(idx_time, datafile + ".time.idx", "T_sorted_index");
		open_index(idx_sys,  datafile + ".sys.idx", "sys_sorted_index");
	}

	void swarmdb::open_index(index_handle &h, const std::string &filename, const std::string &filetype)
	{
		if(access(filename.c_str(), R_OK) != 0)
		{
			ERROR("Cannot open index file '" + filename + "'");
		}

		h.mm.open(filename.c_str(), filetype);
		h.begin = (swarmdb::index_entry *)h.mm.data();
		h.end   = h.begin + h.mm.size()/sizeof(swarmdb::index_entry);
	}
}

//
// Default null-writer (does nothing)
//

class null_writer : public swarm::writer
{
public:
	virtual void process(const char *log_data, size_t length) {}
};

extern "C" swarm::writer *create_writer_null(const std::string &cfg)
{
	return new null_writer();
}

//
// Binary writer
//

class binary_writer : public swarm::writer
{
protected:
	std::auto_ptr<std::ostream> output;
	std::string rawfn, binfn;

public:
	binary_writer(const std::string &cfg);
	virtual void process(const char *log_data, size_t length);
	~binary_writer();
};

extern "C" swarm::writer *create_writer_binary(const std::string &cfg)
{
	return new binary_writer(cfg);
}

binary_writer::binary_writer(const std::string &cfg)
{
	std::string output_fn;
	std::istringstream ss(cfg);
	if(!(ss >> binfn))
		ERROR("Expected 'binary <filename.bin>' form of configuration for writer.")
	rawfn = binfn + ".raw";

	output.reset(new std::ofstream(rawfn.c_str()));
	if(!*output)
		ERROR("Could not open '" + rawfn + "' for writing");

	// write header
	swarm::file_header fh("unsorted_output // Unsorted output file");
	output->write((char*)&fh, sizeof(fh));
}

binary_writer::~binary_writer()
{
	output.reset(NULL);
	if(swarm::sort_binary_log_file(binfn, rawfn))
	{
		unlink(rawfn.c_str());

		// just touch it to auto-generate the indices
		swarm::swarmdb db(binfn);
	}
}

void binary_writer::process(const char *log_data, size_t length)
{
	// TODO: filter out the printfs
	output->write(log_data, length);
}
