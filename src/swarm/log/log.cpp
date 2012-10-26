/*************************************************************************
 * Copyright (C) 2011 by Eric Ford and the Swarm-NG Development Team     *
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

/*! \file log.cpp 
 *    \brief Process the event types. 
 *
 *
 */

#include "log.hpp"
namespace swarm { namespace log {

  int num_ints_for_event(const int code)
  {
    switch(code)
      {
      case EVT_RV_OBS:
      case EVT_ASTROM_OBS:
      case EVT_TIMING_OBS:
      case EVT_DIRECT_IMAGE_OBS:
      case EVT_TRANSIT:
      case EVT_OCCULTATION:
	return 1;
      case EVT_MUTUAL_EVENT:
	return 2;
      default:
	return 0;
      }
  }

  int num_doubles_for_event(const int code)
  {
    switch(code)
      {
      case EVT_RV_OBS:
      case EVT_TIMING_OBS:
	return 1;
      case EVT_ASTROM_OBS:
      case EVT_TRANSIT:
      case EVT_OCCULTATION:
      case EVT_MUTUAL_EVENT:
	return 2;
      case EVT_DIRECT_IMAGE_OBS:
	return 3;
      default:
	return 0;
      }
  }

} } // namespace log:: swarm
