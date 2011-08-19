#pragma once
#include "../swarm_error.h"

namespace swarm {

	struct key_not_found: std::runtime_error {
		key_not_found(const std::string &key): std::runtime_error("Key not found " + key) {}
	};

	class config : public std::map<std::string, std::string> {
	public:
		static config load(const std::string &fn);

		static bool valid_value(const std::string& v, std::string) {
			return v.size() > 0;
		}
		static std::string parse(const std::string& v, std::string) {
			return v;
		}
		static bool valid_value(const std::string& v, double) {
			double d = 0;
			int r = sscanf(v.c_str(),"%lg", &d);
			return r;
		}
		static double parse(const std::string& v, double) {
			double d = 0;
			int r = sscanf(v.c_str(),"%lg", &d);
			return d;
		}
		static bool valid_value(const std::string& v, int) {
			int d = 0;
			int r = sscanf(v.c_str(),"%i", &d);
			return r;
		}
		static int parse(const std::string& v, int) {
			int d = 0;
			int r = sscanf(v.c_str(),"%i", &d);
			return d;
		}

		bool valid(const std::string& key)const{
			return count(key) && (at(key).size() > 0);
		}
		template<typename T> bool valid(const std::string& key , T t )const{
			return valid(key) && valid_value(at(key),t);
		}

		template<typename T> T optional(const std::string& key , const T& default_value)const{
			if( count(key) ) {
				std::string value = at(key);
				return valid_value(value,T()) ? parse(value,T()) : default_value;
			} else
				return default_value;
		}
		template<typename T> T require(const std::string& key, const T&  = T() )const{
			std::string value = at(key);
			if( count(key) && valid_value(value,T()) )
				return parse(value,T());
			else
				throw key_not_found(key);
		}
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


