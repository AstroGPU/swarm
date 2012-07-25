/*************************************************************************
 * Copyright (C) 2008-2010 by Mario Juric & Swarm-NG Development Team    *
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

/*! \file integrator.cpp
 *  \brief class integrator 
 *
*/

#include "common.hpp"
#include "integrator.hpp"
#include "log/logmanager.hpp"
#include "plugin.hpp"
#include "gpu/utilities.hpp"
#include "gpu/device_settings.hpp"

namespace swarm {

	const int integrator::_default_max_iterations = 10000;
	const int integrator::_default_max_attempts   = 1000000;

	void integrator::set_log_manager(log::Pmanager& l){
		_logman = l;
		_log = l->get_hostlog();
	}

        gpulog::host_log* integrator::get_host_log(){
	  //	  return _logman->get_hostlog();
	  return _log;
	}

	void gpu::integrator::set_log_manager(log::Pmanager& l){
		Base::set_log_manager(l);
		set_log(l->get_gpulog());
	}

	integrator::integrator(const config &cfg){
		set_log_manager(log::manager::default_log());
		_max_iterations = cfg.optional("max_iterations", _default_max_iterations );
		_max_attempts = cfg.optional("max_attempts", _default_max_attempts );
	}

	gpu::integrator::integrator(const config &cfg)
		: Base(cfg), _hens(Base::_ens) {
		set_log_manager(log::manager::default_log());
	}

	int number_of_active_systems(defaultEnsemble ens) {
		int count_running = 0;
		for(int i = 0; i < ens.nsys() ; i++)
			if( ens[i].is_active() ) count_running++;
		return count_running;
	}
	int number_of_not_disabled_systems(defaultEnsemble ens) {
		int count_running = 0;
		for(int i = 0; i < ens.nsys() ; i++)
			if( !ens[i].is_disabled() ) count_running++;
		return count_running;
	}
	void activate_all_systems(defaultEnsemble& ens) {
		for(int i = 0; i < ens.nsys() ; i++)
		  {
			if(!ens[i].is_active())
			  ens[i].set_active();
		  }
	}
	void activate_inactive_systems(defaultEnsemble& ens) {
		for(int i = 0; i < ens.nsys() ; i++)
		  {
			if(ens[i].is_inactive())
			  ens[i].set_active();
		  }
	}

	void integrator::integrate() {
		activate_inactive_systems(_ens);
		for(int i = 0; i < _max_attempts; i++)
		  {
			launch_integrator();
			_logman->flush();
			if( number_of_active_systems(_ens) == 0 )
				break;
		}
	};

	void gpu::integrator::integrate() {

		

		activate_inactive_systems(_ens);

		upload_ensemble();
		for(int i = 0; i < _max_attempts; i++)
		  {
			launch_integrator();
			_logman->flush();
			if( number_of_active_systems(_dens) == 0 )
				break;
		}
		download_ensemble();
	};


/*!
   \brief Integrator instantiation support

  @param[in] cfg configuration class
*/

/* Must use factory class to dynamically load integrator subclass
 * instead of using constructor. Done so that users can define their
 * own integrators that the swarm code does not need to know about at
 * compile time
 */
Pintegrator integrator::create(const config &cfg)
{
        Pintegrator integ;

        if(!cfg.count("integrator")) ERROR("Integrator type has to be chosen with 'integrator=xxx' keyword");

        std::string name = cfg.at("integrator");
        std::string plugin_name = "integrator_" + name;

		try {
			integ.reset( (integrator*) plugin::instance(plugin_name,cfg) );
		}catch(plugin_not_found& e){
			ERROR("Integrator " + name + " not found.");
		}

        return integ;
}


}
