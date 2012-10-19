/***************************************************************************
 *   Copyright (C) 2005 by Mario Juric                                     *
 *   mjuric@astro.Princeton.EDU                                            *
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

/*! \file binarystream.hpp
 *   \brief Defines I/O interface in binary form using boost data types. 
 *
 */

/*
 *	This file has been copied/modified from Mario Juric's libpeyton.
 */

#ifndef binarystream_h__
#define binarystream_h__

#define HAVE_BOOST_IOSTREAMS
#ifdef HAVE_BOOST_IOSTREAMS

#include <iostream>
#include <algorithm>
#include <vector>

#include "util.hpp"

#include <boost/type_traits.hpp>

namespace peyton {
namespace io {

	// forward declarations
	template<typename _CharT, typename _Traits = std::char_traits<_CharT> >
		class basic_obstream;
	typedef basic_obstream<char> obstream;
	
	template<typename _CharT, typename _Traits = std::char_traits<_CharT> >
		class basic_ibstream;
	typedef basic_ibstream<char> ibstream;

 	template<typename _CharT, typename _Traits = std::char_traits<_CharT> >
 		class basic_bstream;
 	typedef basic_bstream<char> bstream;

}
}

//
// STL forward declarations
//
namespace std
{
	// Sequences
	template<typename T, typename A> class vector;
	template<typename T, typename A> class deque;
	template<typename T, typename A> class list;

	// Associative Containers
	template<typename T, typename C, typename A> class set;
	template<typename T, typename C, typename A> class multiset;
	template<typename K, typename V, typename C, typename A> class map;
	template<typename K, typename V, typename C, typename A> class multimap;

	// Numerics
	template<typename T> class valarray;
}

//
// Automatically mark a std::pair of PODs as a POD
//
namespace boost
{
	template <typename First, typename Second>
		class is_pod< ::std::pair<First, Second> > : 
			public ::boost::integral_constant<bool, 
				is_pod<First>::value && is_pod<Second>::value
			>
		{};
}

//
// Type traits utilities. Should be separated into it's own header.
//

namespace peyton {
namespace io {

        namespace binary	// auxiliary/debug variables/methods are found here
	{
		// Type information structure for @c manifest
		struct datatype_info
		{
			std::string name;
			int type_size;
			bool ispod;

			datatype_info(const std::string &_name, int _type_size, bool _ispod)
				: name(_name), type_size(_type_size), ispod(_ispod)
				{}
			datatype_info(const datatype_info& di)
				: name(di.name), type_size(di.type_size), ispod(di.ispod)
				{}
		};
		inline bool operator<(const datatype_info& a, const datatype_info& b) { return a.name < b.name; }
		inline bool operator==(const datatype_info& a, const datatype_info& b) { return a.name == b.name; }
		inline std::ostream &operator<<(std::ostream &out, const datatype_info &di);

		typedef std::vector<datatype_info> manifest_t;
		extern manifest_t manifest;	///< A map between @c T and @c sizeof(T)
		extern bool track_manifest; ///< If set to @c true , @c manifest will be populated

		std::ostream &operator<<(std::ostream &out, const manifest_t &manifest);

		/**
			@brief Add typename information to manifest
			
			Adds the type information to @c manifest , but only on
			first invocation for the given type. Used by ibstream::pod_read() and
			obstream::pod_write() methods.
		*/
		template<typename T>
			inline void add_to_manifest()
			{
					static bool passed = false;
					if(track_manifest && !passed)
					{
						std::string tname = peyton::util::type_name<T>();
						typedef typename ::boost::is_pod<T> pod;
						
						// 'T' and 'const T' map to same datatype name, so check if T
						// is already in the vector
						datatype_info dti(tname, sizeof(T), pod::value);
						if(std::find(manifest.begin(), manifest.end(), dti) == manifest.end())
						{
							manifest.push_back(dti);
						}
	
						if(!pod::value)
						{
							std::cerr << "WARNING: Binary I/O as POD for type '" << peyton::util::type_name<T>() << "'\n";
						}
	
						passed = true;
					}
			}
	} // namespace binary

	// Use this to tell stream that your user defined type
	// may be treated as a plain-old-datatype
	#define BLESS_POD(T) \
	        namespace boost {	    \
			template <> \
				class is_pod<T> : public true_type \
					{}; \
		}
	
	
	
	
	//
	// Output stream
	//
	
	template<typename _CharT, typename _Traits>
		class basic_obstream : public std::basic_ostream<_CharT, _Traits>
		{
		protected:
//			explicit basic_obstream() : std::basic_ostream<_CharT, _Traits>() {}
		public:
			template<typename X>
				basic_obstream& write_pod(const X* v, size_t n)
			{
				binary::add_to_manifest<X>();
				this->write(reinterpret_cast<const char *>(v), n*sizeof(X));
				return *this;
			}
		public:
			explicit basic_obstream(std::basic_streambuf<_CharT, _Traits> *sb)
				: std::basic_ostream<_CharT, _Traits>(sb)
				{  }
			explicit basic_obstream(std::basic_ostream<_CharT, _Traits> &out)
				: std::basic_ostream<_CharT, _Traits>(out.rdbuf())
				{  }
		};
	
	//
	// Input stream
	//
	template<typename _CharT, typename _Traits>
		class basic_ibstream : public std::basic_istream<_CharT, _Traits>
		{
		protected:
			explicit basic_ibstream() : std::basic_istream<_CharT, _Traits>() {}
		public:
			template<typename X>
				basic_ibstream& read_pod(X* v, size_t n)
			{
				binary::add_to_manifest<X>();
				this->read((char *)(v), n*sizeof(X));
				return *this;
			}
		public:
			explicit basic_ibstream(std::basic_streambuf<_CharT, _Traits> *sb)
				: std::basic_istream<_CharT, _Traits>(sb)
				{ }
			explicit basic_ibstream(std::basic_istream<_CharT, _Traits> &in)
				: std::basic_istream<_CharT, _Traits>(in.rdbuf())
				{ }
		};

	//
	// Input/output stream
	//
	template<typename _CharT, typename _Traits>
		class basic_bstream : 
			public basic_obstream<_CharT, _Traits>,
			public basic_ibstream<_CharT, _Traits>
		{
		public:
			explicit basic_bstream(std::basic_streambuf<_CharT, _Traits> *sb)
				: basic_obstream<_CharT, _Traits>(sb),
				  basic_ibstream<_CharT, _Traits>(sb)
				{
					//this->init(sb);
				}
			explicit basic_bstream(std::basic_iostream<_CharT, _Traits> &io)
				: basic_obstream<_CharT, _Traits>(io.rdbuf()),
				  basic_ibstream<_CharT, _Traits>(io.rdbuf())
				{
					//this->init(io.rdbuf());
				}
		};

	
	
	
	//
	// Operator declaration macros
	//

#ifdef _MSC_VER
	#define BOSTREAM2(...) peyton::io::obstream& operator<<(peyton::io::obstream &out, __VA_ARGS__)
	#define BISTREAM2(...) peyton::io::ibstream& operator>>(peyton::io::ibstream &in,  __VA_ARGS__)
#else
	#define BOSTREAM2(T...) peyton::io::obstream& operator<<(peyton::io::obstream &out, T)
	#define BISTREAM2(T...) peyton::io::ibstream& operator>>(peyton::io::ibstream &in,  T)
#endif
	
	//
	// Generic POD input/output operators
	//
	
	template<typename T>
		inline BOSTREAM2(const T &v)
		{
			return out.write_pod(&v, 1);
		}
	
	template<typename T>
		inline BISTREAM2(T &v)
		{
			return in.read_pod(&v, 1);
		}
	
	//
	// C and STL string specializations. Note that there's no
	// operator >>(char *), because of possible buffer overflow issues.
	//

	#define RETFAIL(x) if(!(x)) return in;
	
	template<>
	inline BOSTREAM2(char* const& v)
	{
		size_t len = strlen(v);
		out << len;
		return out.write_pod(v, len);
	}
	
	template<>
	inline BOSTREAM2(std::string const& v)
	{
		out << v.size();
		return out.write_pod(v.c_str(), v.size());
	}
	
	template<>
	inline BISTREAM2(std::string &v)
	{
		size_t len;
		RETFAIL(in >> len);
	
		// FIXME: stack overflow likely
		char *buf = (char*)alloca(len+1);
		buf[len] = 0;
		RETFAIL(in.read_pod(buf, len));

		v = buf;
	
		return in;
	}

	namespace details
	{
		template <typename IT>	// unoptimized version
			inline obstream& itwrite(obstream &out, unsigned int size, IT start, const ::boost::false_type&)
			{
				for(IT i = start; size != 0; --size, ++i) { out << *i; }
				return out;
			}
	
		template <typename IT>	// optimized version, for POD arrays
			inline obstream& itwrite(obstream &out, unsigned int size, IT start, const ::boost::true_type&)
			{
				return out.write_pod(start, size);
			}
	}

	template <typename IT>
		inline obstream& itwrite(obstream &out, unsigned int size, IT start)
		{
			out << size;

			// Dispatching to an optimized version, if we're dealing with an
			// array of POD-like types.
			// NOTE: we assume that iterators which are pointers do not have
			// overloaded * and ++ operators. If they do, this code will malfunction.
			typedef ::boost::is_pointer<IT> is_simple;
			typedef ::boost::is_pod<typename std::iterator_traits<IT>::value_type> is_podd;
			typedef ::boost::integral_constant<bool,
				is_simple::value && is_podd::value
				> is_optimizable;

			return details::itwrite(out, size, start, is_optimizable());
		}

	//
	// Reading routines for containers. There are three versions,
	// one optimized for containers linearly stored in memory,
	// the generic one, and one optimized for maps (avoids the
	// unnecessary temporaries of data_type)
	//
	template <typename C>
		inline ibstream& itread(ibstream &in, C &a)
		{
			unsigned int size;
			RETFAIL(in >> size);
	
			a.clear();
	
			typename C::value_type tmp;
			while(size--)
			{
				RETFAIL(in >> tmp);
				a.insert(a.end(), tmp);
			}
			return in;
		}

	template <typename C>
		inline ibstream& itreadvec(ibstream &in, C &a)
		{
			unsigned int size;
			RETFAIL(in >> size);
			a.resize(size);

			typedef ::boost::is_pod<typename C::value_type> is_podd;
			if(is_podd::value)
			{
				in.read_pod(&a[0], size);
			}
			else
			{
				for(int i = 0; i != size; ++i) { in >> a[i]; }
			}
			return in;
		}

	template <typename C>
		inline ibstream& itreadmap(ibstream &in, C &a)
		{
			unsigned int size;
			RETFAIL(in >> size);
	
			a.clear();

			typename C::key_type key;
			while(size--)
			{
				RETFAIL(in >> key);
				RETFAIL(in >> a[key]);
			}
			return in;
		}

	#undef RETFAIL
	
	//
	// STL specializations
	//
	
	template<typename First, typename Second>	// std::pair<First, Second>
		inline BOSTREAM2(const std::pair<First, Second> &v) { return out << v.first << v.second; }
	
	template <typename T, typename A>
		inline BOSTREAM2(const std::vector<T, A> &a) { return itwrite(out, a.size(), &a[0]); }
	template <typename T, typename A>
		inline BOSTREAM2(const std::deque<T, A> &a) { return itwrite(out, a.size(), a.begin()); }
	template <typename T, typename A>
		inline BOSTREAM2(const std::list<T, A> &a) { return itwrite(out, a.size(), a.begin()); }
	
	template <typename T, typename C, typename A>
		inline BOSTREAM2(const std::set<T, C, A> &a) { return itwrite(out, a.size(), a.begin()); }
	template <typename T, typename C, typename A>
		inline BOSTREAM2(const std::multiset<T, C, A> &a) { return itwrite(out, a.size(), a.begin()); }
		
	template <typename K, typename V, typename C, typename A>
		inline BOSTREAM2(const std::map<K, V, C, A> &a) { return itwrite(out, a.size(), a.begin()); }
	template <typename K, typename V, typename C, typename A>
		inline BOSTREAM2(const std::multimap<K, V, C, A> &a) { return itwrite(out, a.size(), a.begin()); }
	
	template <typename T>
		inline BOSTREAM2(const std::valarray<T> &a) { return itwrite(out, a.size(), &a[0]); }
	
	
	
	
	template<typename First, typename Second>	// std::pair<First, Second>
		inline BISTREAM2(std::pair<First, Second> &v) { return in >> v.first >> v.second; }
	
	template <typename T, typename A>
		inline BISTREAM2(std::vector<T, A> &a) { return itreadvec(in, a); }
	template <typename T, typename A>
		inline BISTREAM2(std::deque<T, A> &a) { return itread(in, a); }
	template <typename T, typename A>
		inline BISTREAM2(std::list<T, A> &a) { return itread(in, a); }
	
	template <typename T, typename C, typename A>
		inline BISTREAM2(std::set<T, C, A> &a) { return itread(in, a); }
	template <typename T, typename C, typename A>
		inline BISTREAM2(std::multiset<T, C, A> &a) { return itread(in, a); }
		
	template <typename K, typename V, typename C, typename A>
		inline BISTREAM2(std::map<K, V, C, A> &a) { return itreadmap(in, a); }
	template <typename K, typename V, typename C, typename A>
		inline BISTREAM2(std::multimap<K, V, C, A> &a) { return itreadmap(in, a); }
	
	template <typename T>
		inline BISTREAM2(std::valarray<T> &a) { return itreadvec(in, a); }

} // namespace io
} // namespace peyton

#endif // HAVE_BOOST_IOSTREAMS

using namespace peyton::io;

#endif // binarystream_h__
