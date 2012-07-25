/*************************************************************************
 * Copyright (C) 2011 by Saleh Dindar and the Swarm-NG Development Team  *
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
#include "../common.hpp"

#include "../types/config.hpp"
#include "../plugin.hpp"

#include "io.hpp"
#include "writer.h"

namespace swarm {

template<int NumData>
class event_record
{
public:
  typedef double data_type[NumData];
  int bodid1, bodid2;
  double time;
  data_type data;
  event_record(const int b, const double t, const double* d)
    : bodid1(b), bodid2(-1), time(t)
  {
    for(int i=0;i<NumData;++i) data[i] = d[i];
  }

  event_record(const int b1, const int b2, const double t, const double* d)
    : bodid1(b1), bodid2(b2), time(t)
  {
    for(int i=0;i<NumData;++i) data[i] = d[i];
  }
  int get_num_data() const 
  { return NumData;  };

};

// Specialized version for NumData=-1 to allow for variable length
template<>  class event_record<-1>
{
public:
  typedef std::vector<double> data_type;
  int bodid1, bodid2;
  double time;
  data_type data;
  event_record(const int b, const double t, const data_type& d)
    : bodid1(b), bodid2(-1), time(t), data(d.size())
  {
    for(int i=0;i<data.size();++i) data[i] = d[i];
  }

  event_record(const int b1, const int b2, const double t, const data_type& d)
    : bodid1(b1), bodid2(b2), time(t), data(d.size())
  {
    for(int i=0;i<data.size();++i) data[i] = d[i];
  }

  int get_num_data() const 
  {       return data.size();    }
};

// Specialized version for NumData=0 to allow for no double data (is this needed?)
template<>  class event_record<0>
{
public:
  typedef std::vector<double> data_type;
  int bodid1, bodid2;
  double time;
  //  data_type data;
  event_record(const int b, const double t)
    : bodid1(b), bodid2(-1), time(t)
  {  }

  event_record(const int b1, const int b2, const double t)
    : bodid1(b1), bodid2(b2), time(t)
  {  }

  int get_num_data() const 
  {       return 0;    }
};

/**
 * host_array_writer plugin for use in
 * io.cpp
 *
 */
class host_array_writer : public swarm::writer
{
public: 
  static const int num_doubles_per_event = 2;
  typedef event_record<num_doubles_per_event> event_record_type;
protected:
  std::vector<int> event_codes_to_log;

  typedef std::vector<event_record_type>  event_log_one_system_type;
  typedef std::vector<event_log_one_system_type>  event_log_one_code_type;
  std::vector<event_log_one_code_type> event_log;

public:
  host_array_writer(const config &cfg) : event_codes_to_log(0), event_log(0)
  {
    // eventually figure out how to set multiple event types
    int event_type = cfg.require<int>("host_array_event_type");
    if(event_type!=0) 
      add_event_type_to_log(event_type);
    if(cfg.valid("host_array_event_type2"))
      add_event_type_to_log(cfg.require<int>("host_array_event_type2"));
    if(cfg.valid("host_array_event_type3"))
      add_event_type_to_log(cfg.require<int>("host_array_event_type3"));
    if(cfg.valid("host_array_event_type4"))
      add_event_type_to_log(cfg.require<int>("host_array_event_type4"));
  }
  
  void add_event_type_to_log(const int et)
  {
    event_codes_to_log.push_back(et);
    event_log.resize(event_codes_to_log.size());
  }
  
  event_log_one_code_type& get_event_log_all_systems(const int i)
  {    return event_log[i];  }

  const event_log_one_code_type& get_event_log_all_systems(const int i) const
  {    return event_log[i];  }

  event_log_one_system_type& get_event_log(const int i, const int sys)
  {    return event_log[i][sys];  }

  const event_log_one_system_type& get_event_log_(const int i, const int sys) const
  {    return event_log[i][sys];  }

  int get_event_type(const int i) const
  {    
    if((i>=0) && (i<event_codes_to_log.size()))
      return event_codes_to_log[i];  
    else
      return 0;
  }

  ~host_array_writer()
  {	}
  
  virtual void process(const char *log_data, size_t length)
  {
    int events_added = 0;
    gpulog::ilogstream ils(log_data,length);
    gpulog::logrecord lr;
    while(lr = ils.next()) 
      {
	bool need_to_log = false;
	int event_code = lr.msgid();
	for(int i=0;i<event_codes_to_log.size();++i)
	  {
	    if(event_codes_to_log[i]!=event_code) continue;
	    need_to_log = true;
	    
	    double time;
	    int sysid;
	    lr >> time >> sysid;
	    if(sysid>=event_log[i].size()) event_log[i].resize(sysid+1);
	    
	    // \todo Add code for how to process other observations-type events here
	    //	    if((event_code==EVT_TRANSIT)||(event_code==EVT_OCCULTATION)) 

	    if((event_code>=log::EVT_FIRST_OBS_CODE)&&(event_code<=log::EVT_LAST_OBS_CODE)) 
	      {
		int num_ints = log::num_ints_for_event(event_code);
		int num_doubles = log::num_doubles_for_event(event_code);
		std::vector<int>    intdata(num_ints,-1);
		std::vector<double> doubledata(num_doubles,-1.);
		for(int j=0;j<num_ints;++j)
		  lr >> intdata[j];
		for(int j=0;j<num_doubles;++j)
		  lr >> doubledata[j];
		if(num_doubles>0)
		  {
		    if(num_ints==1)
		      event_log[i][sysid].push_back(event_record_type(intdata[0],time,&doubledata[0]));
		    else if (num_ints==2)
		      event_log[i][sysid].push_back(event_record_type(intdata[0],intdata[1],time,&doubledata[0]));
		  }

		++events_added;
#if DEBUG || 1
		int oldprec = std::cout.precision();
		std::cout.precision(8);
		std::cout << "Event " << event_code << " t= " << time << " sys= " << sysid << " bod= " << intdata[0];
		if(num_ints>=2)
		  std::cout << ", " << intdata[1];
		std::cout.precision(4);
		for(int j=0;j<num_doubles;++j)	  
		  std::cout << " " << doubledata[j];
		std::cout.precision(oldprec);
		std::cout << "\n";
#endif
	      }
	    
	  }
      }
  }
};


writer_plugin_initializer< host_array_writer >
	host_array_writer_plugin("host_array", "This stores selected events in simple arrays on the host");

}

