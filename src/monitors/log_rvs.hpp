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

/*! \file log_rvs.hpp
 *   \brief Defines and implements the monitor that logs time and events at times near a transit on GPU. 
 *
 *
 *  *EXPERIMENTAL*: This class is not thoroughly tested.
 * 
 */

#pragma once

#include<iostream>
#include<sstream>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include "swarm/gpu/gravitation_accjerk.hpp"
#include "swarm/log/gpulog/lprintf.h"

namespace swarm {
  namespace monitors {

#if __CUDA_ARCH__ >= 200

/** Parameters for log_rvs monitor
 * rv_filesname:       name of file with observation times
 * log_rvs_tol (real): anything !=0 turns on logging; otherwise, not really used
 *
 * \ingroup monitors_param
 */ 
struct log_rvs_params {
  double tol, dt;
  int next_time_idx;
  int nbod, num_times;
  thrust::device_ptr<double> rv_times_dev; 

  __host__ log_rvs_params(const config &cfg) 
  {
    num_times = 0; 
    next_time_idx = 0;
    // other paramebers
    tol = cfg.optional("log_rvs_tol", 2.e-8); 
    // WARNING assumes Hermite Fixed w/ nbod bodies
    nbod = cfg.require<int>("nbod"); 
    dt = cfg.require("time_step", 0.0);  // use hermite's timestep

    // Read observation times input device_array
    std::string filename;
    filename = cfg.optional("rv_filename",std::string("rv.in"));
    std::ifstream rv_file(filename.c_str());
    if(!rv_file.good())
      { 
	std::cerr << "# Can't open >" << rv_file << "<.\n"; 
	return ;
      }

    thrust::host_vector<double> rv_times_host;
    while(!rv_file.eof())
      {
	std::string line;
	std::getline(rv_file,line);
	if(rv_file.eof()) continue;
 	std::stringstream line_stream(line);
	double time;
	line_stream >> time;
	try {	rv_times_host.push_back(time); }
	catch(thrust::system_error e)
	  {      std::cerr << "Error push_back: " << e.what() << std::endl;    }
      }
    rv_file.close();
    std::cerr << "# Read " << rv_times_host.size() << " RV observation times from " << filename << "\n";

    try
      { thrust::sort(rv_times_host.begin(),rv_times_host.end()); }
    catch(thrust::system_error e)
      {      std::cerr << "Error sort: " << e.what() << std::endl;    }

    // Allocate & upload list of times onto device
    try
      {	
	num_times = rv_times_host.size();
	rv_times_dev = thrust::device_malloc<double>(num_times);       
      }
    catch(thrust::system_error e)
      { std::cerr << "Error: device_malloc: " << e.what() << std::endl;    }

    try
      { 
/*	for(int i=0;i<rv_times_host.size();++i)
	  {
	    rv_times_dev[i] = rv_times_host[i];
	  } */
		thrust::copy(rv_times_host.begin(),rv_times_host.end(),rv_times_dev);
      }
    catch(thrust::system_error e)
      {      std::cerr << "Error copy: " << e.what() << std::endl;    }


  }


  // WARNING:  This is not automatically called by destructor
  // because I can't figure how to make that happen only on the host destructor
  //
  __host__ void deallocate_device_data()
  {
    thrust::device_free(rv_times_dev);
  }

};

/**
 * Monitor that logs several pieces of information at times near a transit.
 *  *EXPERIMENTAL*: This class is not thoroughly tested.
 *  \ingroup experimental
 * 
 * The information that are recorded are:
 *   * the estimated transit time
 *   * the estimated minimum impact parameter in units of stellar radii 
 *   * the estimated velocity projected onto the planet of the sky in units of stellar radii per time
 * 
 * 
 *   *  Assumes integration results in increasing time.
 *   *  Currently, highly experimental
 *   *  Currently, hardwired for Hermite, Fixed, std. Gravitation, nbod=3
 *   *  
 * 
 *  \ingroup monitors
 *
 *   parameters have yet to be finalized
 */
template<class log_t>
class log_rvs {
	public:
	typedef log_rvs_params params;

	private:
	params _params;
        bool condition_met;

	ensemble::SystemRef& _sys;
  //	double _next_log_time;
	log_t& _log;

	public:
		template<class T>
		static GENERIC int thread_per_system(T compile_time_param){
			return 1;
		}

		template<class T>
		static GENERIC int shmem_per_system(T compile_time_param) {
			 return 0;
		}
        GPUAPI bool is_deactivate_on() { return false; };
        GPUAPI bool is_log_on() { return _params.tol!=0.; };
        GPUAPI bool is_verbose_on() { return false; };
        GPUAPI bool is_any_on() { return is_deactivate_on() || is_log_on() || is_verbose_on() ; }
        GPUAPI bool is_condition_met () { return ( condition_met ); }
        GPUAPI bool need_to_log_system () 
          { return (is_log_on() && is_condition_met() ); }
        GPUAPI bool need_to_deactivate () 
          { return ( is_deactivate_on() && is_condition_met() ); }

        GPUAPI void log_system()  {  log::system(_log, _sys);  }
	
        GPUAPI void operator () (const int thread_in_system) 
          { 
	    if(pass_one(thread_in_system))
	      pass_two(thread_in_system);
	  }

	GPUAPI bool pass_one (const int thread_in_system) 
          {
            condition_met = false;
	    double t = _sys.time();
	    if(_params.next_time_idx>=_params.num_times) return false;
	    while(_params.rv_times_dev[_params.next_time_idx]<t)
	      {
		_params.next_time_idx++;
		if(_params.next_time_idx>=_params.num_times) return false;
	      }
	    double log_time = _params.rv_times_dev[_params.next_time_idx];
	    if( (log_time>=t) && (log_time<t+_params.dt) )
	      return true;
	    else
	      return false;
	  }
	    

  template<int nbod>
  GPUAPI double calc_star_vz(const int& thread_in_system, const double& dt)
  {
    extern __shared__ char shared_mem[];
	typedef swarm::compile_time_params_t<nbod> par_t;
    typedef swarm::gpu::bppt::GravitationAccJerk<par_t> calcForces_t;
    typedef typename calcForces_t::shared_data grav_t;
    calcForces_t calcForces(_sys,*( (grav_t*) (&shared_mem[(swarm::gpu::bppt::sysid_in_block()%SHMEM_CHUNK_SIZE)*sizeof(double)+(swarm::gpu::bppt::sysid_in_block()/SHMEM_CHUNK_SIZE)*SHMEM_CHUNK_SIZE*(nbod*(nbod-1)*3*sizeof(double))]) ) );
    
    double acc, jerk;
#if USING_INTEGRATOR_THAT_OVER_WRITES_SHARED_DATA || 1
    const int bid = swarm::gpu::bppt::thread_body_idx(nbod);
    const int cid = swarm::gpu::bppt::thread_component_idx(nbod);
    double pos = _sys[bid][cid].pos(), vel = _sys[bid][cid].vel();
#if 1
    calcForces(thread_in_system,bid,cid,pos,vel,acc,jerk);
#else
    __syncthreads();
    if(thread_in_system < calcForces_t::pair_count)
     calcForces.calc_pair(thread_in_system);
     __syncthreads();
#endif
#endif
     if(thread_in_system==0)
       {

	 calcForces.sum(0,2,acc,jerk);
	 double star_vz = _sys[0][2].vel()+dt*(acc+0.5*dt*jerk);
	 return star_vz;
       }
     else
       { return 0.; }
  }



	GPUAPI int pass_two (const int thread_in_system) 
          {   
            if(is_any_on())
              {
		double t_log = _params.rv_times_dev[_params.next_time_idx];
		double dt = t_log-_sys.time();
		const int nbod = _params.nbod;
		double star_vz;
		if((nbod==3) && (nbod<MAX_NBODIES))
		  star_vz = calc_star_vz<3>(thread_in_system,dt);
#if 0
		else if((nbod==4) && (nbod<MAX_NBODIES))
		  star_vz = calc_star_vz<4>(thread_in_system,dt);
		else if((nbod==5) && (nbod<MAX_NBODIES))
		  star_vz = calc_star_vz<5>(thread_in_system,dt);
		else if((nbod==6) && (nbod<MAX_NBODIES))
		  star_vz = calc_star_vz<6>(thread_in_system,dt);
		else if((nbod==7) && (nbod<MAX_NBODIES))
		  star_vz = calc_star_vz<7>(thread_in_system,dt);
		else if((nbod==8) && (nbod<MAX_NBODIES))
		  star_vz = calc_star_vz<8>(thread_in_system,dt);
#endif
		if(thread_in_system==0) 
		  {
		    log::event(_log,log::EVT_RV_OBS,t_log,_sys.id(),0,star_vz);	
		  }	
		_params.next_time_idx++;
		condition_met = true;
	      }
	    return condition_met;
	  }


	GPUAPI log_rvs(const params& p,ensemble::SystemRef& s,log_t& l)
	  :_params(p),_sys(s),_log(l) {}
	
};

#endif

  }


}
