#include "common.hpp"
#include "log/io.hpp"
#include <boost/regex.hpp>
#include <boost/program_options.hpp>

namespace swarm { 
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
	
	namespace query {

void execute(const std::string &datafile, swarm::time_range_t T, swarm::sys_range_t sys);
enum planets_coordinate_system_t {
	 astrocentric, barycentric, jacobi
};
void set_keplerian_output(const planets_coordinate_system_t& coordinate_system = astrocentric ) ;

} }

