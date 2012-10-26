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

/*! \file host_array_writer.hpp 
 *    \brief Defines an event writer plug-in for io.cpp.
 *
 *
 *  *EXPERIMENTAL*: This class is not thoroughly tested.
 * 
 *
 */



#include "../common.hpp"
#include "../types/config.hpp"
#include "../plugin.hpp"

#include "io.hpp"
#include "writer.h"

namespace swarm { namespace log {

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

/// Specialized version for NumData=-1 to allow for variable length
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

/// Specialized version for NumData=0 to allow for no double data (is this needed?)
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
 *  A writer plugin that keeps the data in the memory.
 *
 *  *EXPERIMENTAL*: This class is not thoroughly tested.
 * 
 *
 */
class host_array_writer : public writer
{
public: 
  static const int max_num_doubles_per_event = 2;
  typedef event_record<max_num_doubles_per_event> event_record_type;
protected:
  std::vector<int> event_codes_to_log;
  typedef std::vector<event_record_type>  event_log_one_system_type;
  typedef std::vector<event_log_one_system_type>  event_log_one_code_type;
  std::vector<event_log_one_code_type> event_log;
  int debug;

public:
  host_array_writer(const config &cfg);
  
  void add_event_type_to_log(const int et);
  
  event_log_one_code_type& get_event_log_all_systems(const int i)
  {    return event_log[i];  }

  const event_log_one_code_type& get_event_log_all_systems(const int i) const
  {    return event_log[i];  }

  event_log_one_system_type& get_event_log(const int i, const int sys)
  {    return event_log[i][sys];  }

  const event_log_one_system_type& get_event_log(const int i, const int sys) const
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
  
  virtual void process(const char *log_data, size_t length);
};


//writer_plugin_initializer< host_array_writer >
//	host_array_writer_plugin("host_array", "This stores selected events in simple arrays on the host");

} } // namespcae log::swarm

