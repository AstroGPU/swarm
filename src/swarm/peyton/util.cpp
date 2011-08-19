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
#include "util.hpp"

#include <cstring>
#include <cstdio>
#include <cctype>
#include <cstdlib>

#include <string>

using namespace peyton;

char *util::trim(char *dest, const char *src)
{
	int cnt = 0;
	while(isspace(*src) && *src) src++;
	while(*src) dest[cnt++] = *(src++);
	cnt--;
	while(cnt > -1 && isspace(dest[cnt])) cnt--;
	dest[cnt+1] = 0;
	return dest;
}

char *util::trim(char *txt)
{
	char *s = strdup(txt);
	trim(txt, s);
	free(s);
	return txt;
}

///////////

std::string util::ltrim( const std::string &str, const std::string &whitespace)
{

   int idx = str.find_first_not_of(whitespace);
   if( idx != std::string::npos )
       return str.substr(idx);
   
   return "";
}

std::string util::rtrim( const std::string &str, const std::string &whitespace)
{

   int idx = str.find_last_not_of(whitespace);
   if( idx != std::string::npos )
       return str.substr(0,idx+1);

   return str;
}

std::string util::trim( const std::string &str, const std::string &whitespace)
{
    return rtrim(ltrim(str));
}
