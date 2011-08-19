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

#pragma once

#include "types/ensemble.hpp"
#include "types/config.hpp"
#include "log/log.hpp"


namespace swarm {

class integrator;
typedef integrator *(*integratorFactory_t)(const config &cfg);

namespace log {
	class manager;
}

class integrator {	
	protected:
	defaultEnsemble _ens;
	double _destination_time;
	gpulog::host_log* _log;
	const static int _default_max_iterations;
	const static int _default_max_attempts;
	int _max_iterations;
	int _max_attempts;

	public:
	integrator(const config &cfg);
	virtual void launch_integrator() = 0 ;

	virtual void integrate();

	virtual defaultEnsemble& get_ensemble() {
		return _ens;
	}

	virtual void set_ensemble(defaultEnsemble& ens) {
		_ens = ens;
	}

	virtual void set_destination_time(const double& destination_time) {
		_destination_time = destination_time;
	}
	static integrator* create(const config &cfg);

	virtual void set_log_manager(swarm::log::manager& l);

};

namespace gpu {


class integrator : public swarm::integrator {
	typedef swarm::integrator Base;
	protected:
	hostEnsemble& _hens;
	//TODO: use cux auto ptr to make sure we don't have memory leaks
	deviceEnsemble _dens;

	gpulog::device_log* _log;

	public: 

	integrator(const config &cfg);

	virtual void integrate();

	virtual void set_log_manager(log::manager& l);

	void set_log(gpulog::device_log* log) { _log = log; }

	void set_ensemble(defaultEnsemble& ens) {
		_hens = ens;
		_dens = _hens.cloneTo<deviceEnsemble>();
	}

	void upload_ensemble() {
		_hens.copyTo( _dens );
	}
	
	void download_ensemble() {
		_dens.copyTo( _hens );
	}

	virtual dim3 gridDim() = 0;
	virtual dim3 threadDim() = 0;
	virtual int  shmemSize() = 0;
};


}
}
