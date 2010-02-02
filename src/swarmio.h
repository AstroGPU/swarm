#ifndef swarmio_h__
#define swarmio_h__

#include "swarm.h"
#include <astro/binarystream.h>
#include <fstream>
#include <sstream>

// Convert a variable of arbitrary type to a string.
// NOTE: heavy (unoptimized) function, use sparingly
template<typename T>
std::string str(const T& var)
{
	std::ostringstream ss;
	ss << var;
	return ss.str();
}

namespace swarm {

//
// I/O and snapshotting support
//
class ens_writer
{
protected:
	std::string fn;
	std::ofstream out;
	obstream bout;
public:
	ens_writer(const std::string &fn);
	ens_writer &operator <<(const cpu_ensemble &ens);
	operator bool() const { return bout; }
};

class ens_reader
{
protected:
	std::string fn;
	std::ifstream in;
	ibstream bin;
public:
	ens_reader(const std::string &fn);
	ens_reader &operator >>(cpu_ensemble &ens);
	operator bool() const { return bin; }
};

} // end namespace swarm

#endif
