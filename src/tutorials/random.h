/*************************************************************************
 * Copyright (C) 2009-2010 by Eric Ford & the Swarm-NG Development Team  *
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

/*! \file random.h
 *  \brief simple interface to CPU random number generators
 *
 * \todo use better random nubmer generators (e.g., from boost++)
*/

#ifndef H_SWARM_RANDOM

namespace swarm {

double draw_uniform01();
double draw_std_normal();	
double draw_value_from_config(const swarm::config& cfg, const std::string& name, const int bod, double min, double max);

double draw_uniform01() 
{ return static_cast<double>(rand())/static_cast<double>(RAND_MAX); }

double draw_std_normal()
{
  // Box–Muller algorithm
  double rn1 = draw_uniform01();
  double rn2 = draw_uniform01();
  return sqrt(-2.*log(rn1))*cos(2.*M_PI*(rn2));
}

double draw_value_from_config(const swarm::config& cfg, const std::string& name, const int bod, double min, double max)
{
  double val;
  std::stringstream s;

  s.str(""); s << name << '_' << bod << "_min";      
  if(cfg.count(s.str()))
    {
      double min_val = atof(cfg.at(s.str()).c_str());
      min = std::max(min,min_val);
    }
  s.str(""); s << name << '_' << bod << "_max";      
  if(cfg.count(s.str()))
    {
      double max_val = atof(cfg.at(s.str()).c_str());
      max = std::min(max,max_val);
    }

  s.str(""); s << name << '_' << bod;
  if(cfg.count(s.str()))
    {
      val = atof(cfg.at(s.str()).c_str());      
      if((val<min)||(val>max))
	{
	  std::cerr << "# " << s.str() << ": " << val << " ["<< min << ", " << max << "]\n";
	}
      assert(val>=min); 
      assert(val<=max);
      
      s.str(""); s << name << '_' << bod << "_sigma";      
      if(cfg.count(s.str()))
	{
	  double mean = val;
	  double sigma = atof(cfg.at(s.str()).c_str());
	  if(sigma)
	    {
	      do {
		val = mean + sigma*draw_std_normal();
	      } while((val<min)||(val>max));
	    }
	  else
	    { val = mean; }
	}
    }
  else
    {
      val = min + (max-min)*draw_uniform01();
    }
  return val;
}


}

#endif
