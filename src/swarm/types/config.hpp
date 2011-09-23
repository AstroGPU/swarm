#pragma once
#include <stdexcept>

namespace swarm {


	//! Raised when a key not found in the configuration
	struct key_not_found: std::runtime_error {
		key_not_found(const std::string &key): std::runtime_error("Key not found: \"" + key + "\"") {}
	};

	//! Configuration used to configure most classes
	class config : public std::map<std::string, std::string> {
	public:

		config(){}

		/*! Simple constructur for inline creation of config objects form 2D arrays of strings.
		 *  Example:
\code
const char* pairs[][2] = { { "nsys" , "8000" }, { "nbod", "4" } };
config inline_configuration( pairs );
\endcode
		 *
		 *  \args key_value_pairs 2D array of key-value pairs for initializing the config object
		 *
		 */
		template< class T >
		config(const T& key_value_pairs){
			initialize(key_value_pairs, sizeof(T)/sizeof(char*[2]));
		}

		void initialize(const char* const key_value_pairs[][2], const int n);

		/*! Load a new configuration from a text file and merge it with cfg.
		 * Syntax:
\verbatim
key=value
key2=value2
 .
 .
 .
\endverbatim
		 *
		 */
		static config load(const std::string &fn,config cfg = config() );

		//! Is v a valid string
		static bool valid_value(const std::string& v, std::string) {
			return v.size() > 0;
		}
		//! Parse a string
		static std::string parse(const std::string& v, std::string) {
			return v;
		}
		
		//! Is v a valid double
		static bool valid_value(const std::string& v, double) {
			double d = 0;
			int r = sscanf(v.c_str(),"%lg", &d);
			return r == 1;
		}
		//! Parse a double value
		static double parse(const std::string& v, double) {
			double d = 0;
			int r = sscanf(v.c_str(),"%lg", &d);
			return d;
		}

		//! Is v a valid float
		static bool valid_value(const std::string& v, float) {
			float d = 0;
			int r = sscanf(v.c_str(),"%g", &d);
			return r == 1;
		}
		//! Parse a float value
		static float parse(const std::string& v, float) {
			float d = 0;
			int r = sscanf(v.c_str(),"%g", &d);
			return d;
		}

		//! Is v a valid integer
		static bool valid_value(const std::string& v, int) {
			int d = 0;
			int r = sscanf(v.c_str(),"%i", &d);
			return r == 1;
		}
		//! Parse an integer
		static int parse(const std::string& v, int) {
			int d = 0;
			int r = sscanf(v.c_str(),"%i", &d);
			return d;
		}

		//! Query if the value for specified key exists and is a valid entry
		bool valid(const std::string& key)const{
			return count(key) && (at(key).size() > 0);
		}
		//! Query if the value for specified key exists and is a valid value of type T
		//! Example: vaild("log interval", 0.0) makes sure that "log interval" exists and 
		//! is a valid double value ( same type as 0.0 )
		template<typename T> bool valid(const std::string& key , T t )const{
			return valid(key) && valid_value(at(key),t);
		}

		//! Value of an optional parameter key, if it does not exists, then returns the default_value
		template<typename T> T optional(const std::string& key , const T& default_value)const{
			if( count(key) ) {
				std::string value = at(key);
				return valid_value(value,T()) ? parse(value,T()) : default_value;
			} else
				return default_value;
		}
		//! Value of a mandatory parameter key. 
		//! If the key does not exists error is raised. Otherwise text is parsed
		//! as a value of type T and is returned
		template<typename T> T require(const std::string& key, const T&  = T() )const{
			if( count(key) && valid_value(at(key),T()) )
				return parse(at(key),T());
			else
				throw key_not_found(key);
		}
	};

}


