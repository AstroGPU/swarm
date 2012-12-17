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

/*! \file binarystream.cpp
 *   \brief Implements the I/O operators and methods.
 *
 */

 
#define HAVE_BOOST_IOSTREAMS
#ifdef HAVE_BOOST_IOSTREAMS

#include "binarystream.hpp"

#include <iomanip>

std::vector<peyton::io::binary::datatype_info> peyton::io::binary::manifest;
bool peyton::io::binary::track_manifest = true;

std::ostream&
peyton::io::binary::operator<<(std::ostream &out, const datatype_info &di)
{
	return out << std::setw(62) << di.name 
				<< std::setw(8) << di.type_size 
				<< std::setw(8) << di.ispod;
}

std::ostream&
peyton::io::binary::operator<< (std::ostream &out, const manifest_t &manifest)
{
	out << "------------------------------------------------------------------------------\n";
	out << " Binary datatype output manifest                                              \n";
	out << "                                                     data type    size  is_pod\n";
	out << "------------------------------------------------------------------------------\n";
	
	for(int i = 0; i != manifest.size(); i++)
	{
		out << manifest[i] << "\n";
	}
	return out;
}

#include <fstream>
#include <map>
#include <valarray>

#pragma pack(1)
struct test_t {
	char c;
	short s;
	float f;
	double d;
};
#pragma pack()
BLESS_POD(test_t);

int demo_binarystream()
{
	using namespace peyton::io;
	using namespace std;

	// (Not quite) typical usage with an iostream, both reading and writing
	// A more typical usage would be with istream/ibstream (for input) and 
	// ostream/obstream (for output). Alternatively, the ?bstream classes'
	// constructor accepts std::streambufs, so you can use things like
	// boost::iostreams to attain filtering, compression, etc...
	fstream io("bla.bin", ios::in | ios::out | ios::binary | ios::trunc);
	bstream bio(io);
	bio.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit);

	// This writes a valarray to disk
	valarray<int> v, w;
	v.resize(20);
	FOR(0, v.size()) { v[i] = i; }
	bio << v;
	// This appends a plain structure
	test_t t = { 'a', -2, 1., 2. }, x;
	bio << t;
	// This serializes a std::map<int, string>
	typedef map<int, std::string> M;
	M a, b;
	a[1] = "jedan"; a[2] = "dva"; a[3] = "tri";
	bio << a;

	// This reads the valarray back
	bio.seekg(0);
	bio >> w;
	cout << "number of elements read: " << w.size() << "\nelements:";
	FOR(0, w.size()) { cout << " " << w[i]; }
	cout << "\n\n";
	// This reads the POD-like structure back
	bio >> x;
	cout << "User defined plain struct: " << x.c << " " << x.s << " " << x.f << " " << x.d << "\n\n";
	// This reads the map back in
	bio >> b;
	cout << "map size: " << b.size() << "\n";
	cout << "map elements:";
	FOREACH2(M::const_iterator, b) { cout << " (" << (*i).first << "->" << (*i).second << ")"; }
	cout << "\n\n";

	// This prints out all datatypes used in binary output, together
	// with their size (given by sizeof()), and weather have they been
	// declared as PODs (plain-ol'-datatypes) or not. Good for debugging.
	// Collection of this information can be turned off by setting
	// peyton::io::binary::track_manifest to false. It's turned on by
	// default.
	cout << binary::manifest << "\n";

	return 0;
}

#if 0
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

void demo_binarystream_with_boost_iostreams()
{
	using namespace boost::iostreams;
	using namespace peyton::io;

	filtering_streambuf<output> out(bzip2_compressor() | file_sink("bla.bin.bz2"));
	obstream bout(&out);

	std::valarray<int> v, w;
	v.resize(5);
	FOR(0, v.size()) { v[i] = i; }
	bout << v;
	out.reset();

	filtering_streambuf<input> in(bzip2_decompressor() | file_source("bla.bin.bz2"));
	ibstream bin(&in);
	bin >> w;
	std::cout << "size=" << w.size() << "\n";
	FOR(0, w.size()) { std::cout << w[i] << " "; } std::cout << "\n";
}
#endif

#else // HAVE_BOOST_IOSTREAMS

#include <iostream>

int demo_binarystream()
{
	std::cerr << "Boost.Iostream support not compiled in.";
	return -1;
}

#endif // HAVE_BOOST_IOSTREAMS
