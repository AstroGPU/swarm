/***************************************************************************
 *   Copyright (C) 2005 by Mario Juric   *
 *   mjuric@astro.Princeton.EDU   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*! \file util.hpp
 *   \brief Defines string operation utility functions. 
 *
 */

#ifndef _astro_util_h
#define _astro_util_h

#include "macros.h"
#include "types.hpp"
#include "constants.hpp"

#include <cmath>
#include <string>
#include <cstdlib>
#include <memory.h>
#include <ctype.h>
#include <stdio.h>
#include <exception>
#include <stdexcept>

#ifndef _WIN32
#include <cxxabi.h>
#endif

#ifndef NULL
#define NULL 0
#endif

namespace peyton {
// contains string processing functions, type name demangling,
namespace util {

	/// return approximate longitude of the Sun for a given \a time
	//	Radians approxSunLongitude(peyton::MJD t);

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
	inline std::string str(size_t n) { char buf[20]; sprintf(buf, "%lu", n); return buf; }
	/// convert int to std::string
	inline std::string str(int n) { char buf[20]; sprintf(buf, "%d", n); return buf; }
	/// convert char to std::string
	inline std::string str(char c) { char buf[2] = {c, 0}; return buf; }
	/// convert double to std::string
	inline std::string str(double n, const char *fmt = "%f") { char buf[20]; sprintf(buf, fmt, n); return buf; }
	#if 0
	/// Convert a variable of arbitrary type to a string.
	/// NOTE: heavy (unoptimized) function, use sparringly
	template<typename T>
	inline std::string str(const T& var)
	{
		std::ostringstream ss;
		ss << var;
		return ss.str();
	}
	#endif

	struct demangling_error : public std::runtime_error {
		demangling_error(const char* s): std::runtime_error(("While demangling " + std::string(s)).c_str()) {}
	};

#ifdef _WIN32
	inline std::string type_name(const std::type_info &ti){
		// TODO: find out how to do demangling in Windows
		return ti.name();
	}
#else

	inline std::string type_name(const std::type_info &ti){
		int status = 0;
		char * tmp = abi::__cxa_demangle(ti.name(), 0, 0, &status);
		std::string name = tmp;
		free(tmp);
		if(status != 0)
			throw demangling_error(ti.name());
		return name;
	}

#endif

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

using namespace peyton::util;

#endif
