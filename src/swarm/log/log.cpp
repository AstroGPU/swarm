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
