#include <astro/util.h>

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
