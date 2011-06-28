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

#include <limits>
#include "hp/ensemble.hpp"

namespace swarm {
namespace hp {


class stop_on_ejection {
	double _rmax_squared;

	public: 
	class tester {
		const double & _rmax_squared;
		ensemble::SystemRef& _sys;

		public:

		GPUAPI tester(const float& rmax_squared, ensemble::SystemRef& sys)
			:_rmax_squared(rmax_squared),_sys(sys){}

		GPUAPI bool operator () () { 
			for(int b = 1 ; b < _sys.nbod(); b ++ ){
				if(_sys.distance_squared_between(b,0) > _rmax_squared )
					return true;
			}
			return false; 
		}
	};

	stop_on_ejection(const config &cfg)
	{
		if(!cfg.count("rmax"))
			_rmax_squared = std::numeric_limits<float>::max();
		else
			_rmax_squared = sqr (atof(cfg.at("rmax").c_str()) ); 
	}

	GPUAPI tester get_tester(ensemble::SystemRef& sys){
		return tester(_rmax_squared,sys);
	}
	
};

}
}
