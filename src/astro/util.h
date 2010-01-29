#ifndef _astro_util_h
#define _astro_util_h

#include <astro/macros.h>
#include <astro/types.h>
#include <astro/constants.h>
#include <cmath>
#include <string>
#include <cstdlib>
#include <memory.h>
#include <ctype.h>
#include <stdio.h>
#include <cxxabi.h>

#ifndef NULL
#define NULL 0
#endif

namespace peyton {
/// everything else
namespace util {

	/// return approximate longitude of the Sun for a given \a time
	Radians approxSunLongitude(peyton::MJD t);

	/// remove whitespace from the begining of the string \a str
	std::string ltrim(const std::string &str, const std::string &whitespace = "\t ");
	/// remove whitespace from the end of the string \a str
	std::string rtrim(const std::string &str, const std::string &whitespace = "\t ");
	/// remove whitespace from the begining and the end of the string \a str
	std::string  trim(const std::string &str, const std::string &whitespace = "\t ");
	/** \deprecated Use the version which takes std::string as an argument */
	char *trim(char *txt);
	/** \deprecated Use the version which takes std::string as an argument */
	char *trim(char *dest, const char *src);

	/// pad to given number of characters
	std::string pad(const std::string &s, size_t n, char c = ' ');

	/// Convert all occurences of \\" and \\' in a string to " and '
	std::string unescape(const std::string &str);

	/// convert string to lowercase
	inline std::string tolower(const std::string &s) { std::string o(s); FOREACH2(std::string::iterator, o) { *i = ::tolower(*i); }; return o; }
	/// convert string to uppercase
	inline std::string toupper(const std::string &s) { std::string o(s); FOREACH2(std::string::iterator, o) { *i = ::toupper(*i); }; return o; }

	/// convert size_t to std::string
	inline std::string str(size_t n) { char buf[20]; sprintf(buf, "%d", n); return buf; }
	/// convert int to std::string
	inline std::string str(int n) { char buf[20]; sprintf(buf, "%d", n); return buf; }
	/// convert char to std::string
	inline std::string str(char c) { char buf[2] = {c, 0}; return buf; }
	/// convert double to std::string
	inline std::string str(double n, const char *fmt = "%f") { char buf[20]; sprintf(buf, fmt, n); return buf; }

	/// type-name demangler
	inline std::string type_name(const std::type_info &ti)
	{
		std::string tmp;
		int status;
	
		char *name = abi::__cxa_demangle(ti.name(), 0, 0, &status);
		if(status != 0)
		{
			// TODO: thrown an exception here? If so, make sure to free(name) first!
			tmp = "demangling error (status code = ";
			tmp += str(status);
		}
		else
		{
			tmp = name;
		}
		free(name);

		return tmp;
	}

	/// type-name demangler
	template<typename T>
		std::string type_name()
		{
			return type_name(typeid(T));
		}

	/// type-name demangler
	template<typename T>
		std::string type_name(const T& t)
		{
			return type_name<T>();
		}
}
namespace Util = util; // backwards compatibility hack

}

#define __peyton_util peyton::util

#endif
