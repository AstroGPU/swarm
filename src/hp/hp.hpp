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

#include <cuda_runtime_api.h>
#include "swarm.h"
#include "datatypes.hpp"
#include "ensemble.hpp"


namespace swarm {
namespace hp {

class integrator {	
	protected:
	defaultEnsemble _ens;
	double _destination_time;

	public:
	integrator(const config &cfg){}
	virtual void integrate() = 0 ;

	virtual defaultEnsemble& get_ensemble() {
		return _ens;
	}

	virtual void set_ensemble(defaultEnsemble& ens) {
		_ens = ens;
	}

	virtual void set_duration(const double& duration) {
		_destination_time = duration;
	}

};

	defaultEnsemble convert_ensemble(cpu_ensemble& ens){
		hostEnsemble _hens = hostEnsemble::create(ens.nbod(),ens.nsys());

		// Copy all
		for(int i= 0; i < ens.nsys() ; i++){
			for(int j = 0; j < ens.nbod(); j++) {
				_hens[i][j][0].pos() = ens[i][j].p(0);
				_hens[i][j][1].pos() = ens[i][j].p(1);
				_hens[i][j][2].pos() = ens[i][j].p(2);
				_hens[i][j][0].vel() = ens[i][j].v(0);
				_hens[i][j][1].vel() = ens[i][j].v(1);
				_hens[i][j][2].vel() = ens[i][j].v(2);
				_hens[i][j].mass() = ens[i][j].mass();
			}
			_hens[i].time() = ens[i].time();
		}

		return _hens;
	}

namespace gpu {


class integrator : public hp::integrator {
	typedef hp::integrator Base;
	protected:
	hostEnsemble& _hens;
	//TODO: use cux auto ptr to make sure we don't have memory leaks
	deviceEnsemble _dens;

	public: 

	integrator(const config &cfg): Base(cfg), _hens(Base::_ens) {}

	void integrate() {
		upload_ensemble();
		launch_integrator();
		download_ensemble();
	}

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

	virtual void launch_integrator() = 0;

	virtual dim3 gridDim() = 0;
	virtual dim3 threadDim() = 0;
	virtual int  shmemSize() = 0;
};

}
}
}
