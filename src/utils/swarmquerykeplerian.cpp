/*************************************************************************
 * Copyright (C) 2011 by Mario Juric and the Swarm-NG Development Team   *
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

/*! \file swarmquerykeplerian.cpp
 *   \brief Implements a program for extracting text information from binary files generated 
 *          by swarm GPU logging subsystem. 
 *
 * swarmquery is a useful tool for extracting information about selected systems at all or selected times.
 * usage: swarmquery datafile [Options]
 * swarmquery --help for overview of options:
 *  --help                produce help message
 * -t [ --time ] arg     range of times to query
 * -s [ --system ] arg   range of systems to query
 * --datafile arg        the datafile to query
 * 
 * A range has the format xx..yy.  
 * Since system id's are integers, the range 
 * for the system id, can also be a single integer.
 * 
 * Example:
 * ../bin/swarmquery -s 42 -t 0.1..0.2002 log.bin 
 *
*/

#define ASTROCENTRIC 1
#define BARRYCENTRIC 0 // Not implemented
#define JACOBI 0       // Not implemented

#include <iostream>
#include <string>
#include <sstream>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>

#include "swarm/io.hpp"
#include "tutorials/kepler.h"

// using namespace boost;
// using namespace boost::program_options;

namespace swarm
{
	template<typename T>
	T arg_parse(const std::string &s)
	{
		if(s == "MIN") { return MIN; }
		if(s == "MAX") { return MAX; }
		return boost::lexical_cast<T>(s);
	}

	//
	// Parse ranges of the form:
	//	<r1>..<r2>
	//	MIN..<r2>
	//	<r1>..MAX
	//	ALL
	//
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

	template<typename T>
	std::ostream &operator<<(std::ostream &out, const swarm::range<T> &r)
	{
		if(r.first == T(MIN) && r.last == T(MAX)) { return out << "ALL"; }
		if(r.first == r.last) { return out << r.first; }

		if(r.first == T(MIN)) { out << "MIN"; } else { out << r.first; }
		out << "..";
		if(r.last == T(MAX)) { out << "MAX"; } else { out << r.last; }

		return out;
	}
}

void get_Tsys(gpulog::logrecord &lr, double &T, int &sys)
{
	//std::cerr << "msgid=" << lr.msgid() << "\n";
	if(lr.msgid() < 0)
	{
		// system-defined events that have no (T,sys) heading
		T = -1; sys = -1;
	}
	else
	{
		lr >> T >> sys;
	}
}

// Default output, if no handler is registered
std::ostream& record_output_default(std::ostream &out, gpulog::logrecord &lr)
{
	double T; int sys;
	get_Tsys(lr, T, sys);
	out << lr.msgid() << "\t" << T << "\t" << sys;
	return out;
}

// EVT_SNAPSHOT
std::ostream& record_output_1(std::ostream &out, gpulog::logrecord &lr)
{
	double T;
	int nbod, sys, flags;
	const swarm::body *bodies;
	lr >> T >> sys >> flags >> nbod >> bodies;

	size_t bufsize = 1000;
	char buf[bufsize];
	// arrays with Keplerian orbital elements of body's >=1
	std::vector<double> a(nbod), e(nbod), i(nbod), w(nbod), O(nbod), M(nbod);

	// coordinates of center for each orbit (changes per body in Jacobi)
#if ASTROCENTRIC
	const swarm::body &bstar = bodies[0];
	double cx = bstar.x;
	double cy = bstar.y;
	double cz = bstar.z;
	double cvx = bstar.vx;
	double cvy = bstar.vy;
	double cvz = bstar.vz;
	double cmass = bstar.mass;
#else
#if BARRYCENTRIC
	double cx = 0., cy=0., cz=0., cvx=0., cvy=0., cvz=0., cmass=0.;
	for(int bod = 0; bod != nbod; bod++)
	{
		const swarm::body &b = bodies[bod];
	  cx += b.mass*b.x;
	  cy += b.mass*b.y;
	  cz += b.mass*b.z;
	  cvx += b.mass*b.vx;
	  cvy += b.mass*b.vy;
	  cvz += b.mass*b.vz;
	  cmass += b.mass;
	}
	cx  /= cmass;
	cy  /= cmass;
	cz  /= cmass;
	cvx /= cmass;
	cvy /= cmass;
	cvz /= cmass;
#else
#if JACOBI
	double cx = 0., cy=0., cz=0., cvx=0., cvy=0., cvz=0.;
	double cmass=bstar.mass;
#else
	#error "You have to use some coordinate system!";
#endif
#endif
#endif
	
	for(int bod = 1; bod != nbod; bod++)  // Skip star since printing orbits
	{
	  double x, y, z, vx, vy, vz, mass;
	  const swarm::body &b = bodies[bod];
#if JACOBI
	  cx += b.mass*b.x;
	  cy += b.mass*b.y;
	  cz += b.mass*b.z;
	  cvx += b.mass*b.vx;
	  cvy += b.mass*b.vy;
	  cvz += b.mass*b.vz;
	  cmass += b.mass;

	  y = b.x-cx/cmass;
	  y = b.y-cy/cmass;
	  z = b.z-cz/cmass;
	  vx = b.vx-cx/cmass;
	  vy = b.vy-cy/cmass;
	  vz = b.vz-cz/cmass;
#else
	  x = b.x-cx;
	  y = b.y-cy;
	  z = b.z-cz;
	  vx = b.vx-cvx;
	  vy = b.vy-cvy;
	  vz = b.vz-cvz;
	  cmass = bstar.mass + b.mass;
#endif
	  calc_keplerian_for_cartesian(a[bod], e[bod], i[bod], O[bod], w[bod], M[bod], x, y, z, vx, vy, vz, cmass);

	  const double rad2deg = 180./M_PI;

	  if(bod != 0) { out << "\n"; }
	  snprintf(buf, bufsize, "%10d %g  %5d %5d  %g  % 9.5f % 9.5f % 9.5f  % 9.5f % 9.5f % 9.5f  %d", lr.msgid(), T, sys, bod, b.mass, a[bod], e[bod], i[bod]*rad2deg, O[bod]*rad2deg, w[bod]*rad2deg, M[bod]*rad2deg, flags);
	  out << buf;
	}
	return out;
}

// EVT_EJECTION
std::ostream& record_output_2(std::ostream &out, gpulog::logrecord &lr)
{
	double T;
	int sys;
	swarm::body b;
	lr >> T >> sys >> b;

        size_t bufsize = 1000;
        char buf[bufsize];
     snprintf(buf, bufsize, "%10d %g  %5d %5ld  %g  % 9.5f % 9.5f % 9.5f  % 9.5f % 9.5f % 9.5f", lr.msgid(), T, sys, b.bod, b.mass, b.x, b.y, b.z, b.vx, b.vy, b.vz);
	out << buf;

	return out;
}

typedef std::ostream& (*record_output_function_t)(std::ostream &out, gpulog::logrecord &lr);
record_output_function_t record_output_function_pointers[] = { 0, record_output_1, record_output_2 };

std::ostream &output_record(std::ostream &out, gpulog::logrecord &lr)
{
	int evtid = lr.msgid();

	record_output_function_t fun = record_output_function_pointers[evtid];
	fun = fun ? fun : record_output_default;

	return fun(out, lr);
}

void execute_query(const std::string &datafile, swarm::time_range_t T, swarm::sys_range_t sys)
{
	using namespace swarm;

	swarmdb db(datafile);
	swarmdb::result r = db.query(sys, T);
	gpulog::logrecord lr;
	while(lr = r.next())
	{
		output_record(std::cout, lr);
		std::cout << "\n";
	}
}

int main(int argc, char **argv)
{
	namespace po = boost::program_options;
	using namespace swarm;

	po::options_description desc(std::string("Usage: ") + argv[0] + " <datafile>\n\nOptions");
	desc.add_options()
		("help", "produce help message")
		("time,t", po::value<time_range_t>(), "range of times to query")
		("system,s", po::value<sys_range_t>(), "range of systems to query")
		("datafile", po::value<std::string>(), "the datafile to query")
	;
	po::positional_options_description pd;
	pd.add("datafile", 1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().positional(pd).run(), vm);
	po::notify(vm);

	if (vm.count("help") || !vm.count("datafile")) { std::cout << desc << "\n"; return 1; }

	time_range_t T;
	sys_range_t sys;
	if (vm.count("time")) { T = vm["time"].as<time_range_t>(); }
	if (vm.count("system")) { sys = vm["system"].as<sys_range_t>(); }

        std::cout.flush(); 
	std::cerr << "Printing outputs satisfying T=" << T << " sys=" << sys << "\n";

	std::string datafile(vm["datafile"].as<std::string>());

	execute_query(datafile, T, sys);
}



