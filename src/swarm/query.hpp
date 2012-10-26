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

/*! \file query.hpp
 *   \brief Defines class arg_parse and validate function for data values.
 *
 *
*/


#include "common.hpp"
#include "log/io.hpp"
#include <boost/regex.hpp>
#include <boost/program_options.hpp>

namespace swarm { namespace query {
	
template<typename T>
T arg_parse(const std::string &s)
{
	if(s == "MIN") { return MIN; }
	if(s == "MAX") { return MAX; }
	return boost::lexical_cast<T>(s);
}

/*!
 * Parser for range datatype to use with boost::program_options.
 * 
 * Parse ranges of the form:
 *	<r1>..<r2>
 *	MIN..<r2>
 *	<r1>..MAX
 *	ALL
 * 
 */
template<typename T>
void validate(boost::any& v,
	const std::vector<std::string>& values,
	range<T>* target_type, int)
{
	using namespace boost::program_options;
	typedef range<T> rangeT;

	// Make sure no previous assignment to 'a' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	const std::string& s = validators::get_single_string(values);

	if(s == "ALL")
	{
		v = boost::any(rangeT(ALL));
		return;
	}

	static boost::regex r("(.+?)(?:\\.\\.(.+))?");
	boost::smatch match;
	if (boost::regex_match(s, match, r)) {
		//for(int i=0; i != match.size(); i++) { std::cerr << match[i] << "\n"; }
		if(match[2] == "")
		{
			v = boost::any(rangeT(arg_parse<T>(match[1])));
		}
		else
		{
			v = boost::any(rangeT(arg_parse<T>(match[1]), arg_parse<T>(match[2])));
		}
	} else {
#if BOOST_VERSION  < 104200
		throw validation_error("invalid value");
#else
		throw validation_error(validation_error::invalid_option_value);
#endif
	}
}

/**
 * Petty print a range object
 * 
 */
template<typename T>
std::ostream &operator<<(std::ostream &out, const range<T> &r)
{
	if(r.first == T(MIN) && r.last == T(MAX)) { return out << "ALL"; }
	if(r.first == r.last) { return out << r.first; }

	if(r.first == T(MIN)) { out << "MIN"; } else { out << r.first; }
	out << "..";
	if(r.last == T(MAX)) { out << "MAX"; } else { out << r.last; }

	return out;
}
	

/*! Execute a query on the datafile. The query consists of ranges for time (T), systems (sys) and bodies (bod)
 * 
 * @param datafile Filename for the swarm binary log to query from
 * @param T        Range of times to output
 * @param sys      Range of systems to output
 * @param bod      Range of bodies to output
 * 
 */
void execute(const std::string &datafile, time_range_t T, sys_range_t sys, body_range_t bod = body_range_t() );

enum planets_coordinate_system_t {
  astrocentric, barycentric, jacobi, origin
};

//! Set the output format of the @ref execute to be in Cartesian coordinates
void set_cartesian_output(const planets_coordinate_system_t& coordinate_system = origin) ;
//! Set the output format of the @ref execute to be in Keplerian coordinates
void set_keplerian_output(const planets_coordinate_system_t& coordinate_system = jacobi) ;
/*! Set the output format of the @ref execute for the coordinates. 
 * @param coordinate_system Choice of coordinate system, c.f. @ref planets_coordinate_system_t
 * 
 */
void set_coordinate_system(const planets_coordinate_system_t& coordinate_system);  

} }

