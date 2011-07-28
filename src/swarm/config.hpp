#pragma once
#include "swarm/swarm_error.h"
#include <map>
#include <string>
namespace swarm {

	class config : public std::map<std::string, std::string> {
	};

	/*!
	  \brief get a configuration value for 'key', throwing an error if it doesn't exist

NOTE: heavy (unoptimized) function, use sparingly
*/
	template<typename T>
		void get_config(T &val, const config &cfg, const std::string &key)
		{
			if(!cfg.count(key)) { ERROR("Configuration key '" + key + "' missing."); }
#ifndef __CUDACC__ // CUDA 2.2 C++ bug workaround
			std::istringstream ss(cfg.at(key));
			ss >> val;
#endif
		}

}


