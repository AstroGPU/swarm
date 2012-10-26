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

/*! \file host_array_writer.cpp 
 *    \brief Processes the event log data.
 *
 *
 */

#include "host_array_writer.hpp"

#include "../common.hpp"

#include "../types/config.hpp"
#include "../plugin.hpp"

//#include<unordered_map>
#include "io.hpp"
#include "writer.h"

namespace swarm { namespace log {

  host_array_writer::host_array_writer(const config &cfg) : event_codes_to_log(0), event_log(0)
  {
    debug = cfg.optional<int>("debug_host_array_writer",0);
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
  
  void host_array_writer::add_event_type_to_log(const int et)
  {
    //    int size = event_codes_to_log.size();
    //    event_code_index[et] = size;
    event_codes_to_log.push_back(et);
    event_log.resize(event_codes_to_log.size());
  }
  
  void host_array_writer::process(const char *log_data, size_t length)
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

	    //	 int event_id = (depth>0.) ? log::EVT_TRANSIT : log::EVT_OCCULTATION;
	    //	 log::event(_log,event_id,_sys.time()+dt_min_b2,_sys.id(),j,b,vproj);

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
		if(debug)
		  {
		    int oldprec = std::cout.precision();
		    std::cout.precision(10);
		    std::cout << "EventCode= " << event_code << " t= " << time << " sys= " << sysid << " bod= " << intdata[0];
		    if(num_ints>=2)
		      {
			for(int j=1;j<num_ints;++j)
			  std::cout << ", " << intdata[j];
		      }
		    std::cout.precision(6);
		    for(int j=0;j<num_doubles;++j)	  
		      std::cout << " " << doubledata[j];
		    std::cout.precision(oldprec);
		    std::cout << " \n";
		  }
	      }
	    
	  }
      }
  }


writer_plugin_initializer< host_array_writer >
	host_array_writer_plugin("host_array", "This stores selected events in simple arrays on the host");

} } // namespcae log::swarm

