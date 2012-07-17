//
//
// Author: Mario Juric <mjuric@cfa.harvard.edu>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//

/*! \file query.cpp
 *    \brief routines for extracting text information from binary files generated by swarm GPU logging subsystem.
 * @authors Mario Juric
 *
 * 
 * A range has the format xx..yy.  
 * Since system id's are integers, the range 
 * for the system id, can also be a single integer.
 * 
 *
*/

#include "query.hpp"
#include "kepler.h"


namespace swarm { namespace query {



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

bool keplerian_output = false;
planets_coordinate_system_t planets_coordinate_system = jacobi;

void set_keplerian_output(const planets_coordinate_system_t& coordinate_system) 
{  keplerian_output = true;  planets_coordinate_system = coordinate_system; }

void set_cartesian_output(const planets_coordinate_system_t& coordinate_system)
{  keplerian_output = false;  planets_coordinate_system = coordinate_system; }

void set_coordinate_system(const planets_coordinate_system_t& coordinate_system)  
{  planets_coordinate_system = coordinate_system; }


struct keplerian_t {
	double a, e, i, O, w , M;
};

keplerian_t keplerian_for_cartesian ( const swarm::body& b, const swarm::body& center ) 
{ 
  keplerian_t k;
  double x = b.x - center.x;
  double y = b.y - center.y;
  double z = b.z - center.z;
  double vx = b.vx - center.vx;
  double vy = b.vy - center.vy;
  double vz = b.vz - center.vz;
  double mass = b.mass + center.mass;
  calc_keplerian_for_cartesian(k.a, k.e, k.i, k.O, k.w, k.M, x, y,z, vx, vy, vz, mass );
  return k;
}

body center_of_mass(const body* bodies, const int nbod ){ 
	body center;
	center.x = center.y = center.z = center.vx = center.vy = center.vz = 0.;
	center.mass = 0.;
	int i=0;
	while(i<nbod)
	  {
	    center.x += bodies[i].x*bodies[i].mass;
	    center.y += bodies[i].y*bodies[i].mass;
	    center.z += bodies[i].z*bodies[i].mass;
	    center.vx += bodies[i].vx*bodies[i].mass;
	    center.vy += bodies[i].vy*bodies[i].mass;
	    center.vz += bodies[i].vz*bodies[i].mass;
	    center.mass += bodies[i].mass;
	    i++;
	  }
	center.x /=   center.mass;
	center.y /=   center.mass;
	center.z /=   center.mass;
	center.vx /=  center.mass;
	center.vy /=  center.mass;
	center.vz /= center.mass;

	return center;
}

// EVT_SNAPSHOT
    std::ostream& record_output_1(std::ostream &out, gpulog::logrecord &lr, swarm::body_range_t &body_range)
{
	double time;
	int nbod, sys, flags;
	const body *bodies;
	lr >> time >> sys >> flags >> nbod >> bodies;

	body center;
	const swarm::body &star = bodies[0];

	    switch(planets_coordinate_system)
	      {
	      case astrocentric:
		center = star;
		break;
	      case barycentric:
		center = center_of_mass( bodies, nbod );
		break;
	      case jacobi:
		center = star;
		break;
	      case origin:
		center.x = center.y = center.z = center.vx = center.vy = center.vz = 0.; center.mass = star.mass;
		break;
	      };

	if((time<=0.) && false)  // was used for debugging at some point
	  {
	    if (keplerian_output)
	      { std::cerr << "# Output in Keplerian coordinates  "; }
	    else { std::cerr << "# Output in Cartesian coordinates  "; }
	    switch(planets_coordinate_system)
	      {
	      case astrocentric:
		std::cerr << "(astrocentric) " << center.x << ' ' << center.y << ' '<< center.z << "  " << center.vx << ' ' << center.vy << ' ' << center.vz;
		break;
	      case barycentric:
		std::cerr << "(barycentric) "<< center.x << ' ' << center.y << ' '<< center.z << "  " << center.vx << ' ' << center.vy << ' ' << center.vz;
		break;
	      case jacobi:
		std::cerr << "(jacobi) " << center.x << ' ' << center.y << ' '<< center.z << "  " << center.vx << ' ' << center.vy << ' ' << center.vz;
		break;
	      case origin:
		std::cerr << "(origin) " << center.x << ' ' << center.y << ' '<< center.z << "  " << center.vx << ' ' << center.vy << ' ' << center.vz;
		break;
	      }
	    std::cerr << "\n";
	  }

	
	size_t bufsize = 1000;
	char buf[bufsize];
	for(int bod = 0; bod < nbod; bod++)
	{
	  if(!body_range.in(bod)) { continue; }
	  const swarm::body &b = bodies[bod];
	  if(  keplerian_output && (bod==0) ) { continue; }
	  if(  (planets_coordinate_system==jacobi) && (bod==0) ) { continue; }
	  if( ( !keplerian_output && (planets_coordinate_system!=jacobi) && (bod> 0) ) ||  
	      ( !keplerian_output && (planets_coordinate_system==jacobi) && (bod> 1) ) ||  
	      (  keplerian_output && (bod> 1 ) ) ){ out << "\n"; }
	  
	  if( planets_coordinate_system == jacobi )
	    { center = center_of_mass ( bodies, bod); }
	  
	  if( keplerian_output  && bod > 0) 
	    {
	      if(planets_coordinate_system==barycentric) center.mass -= b.mass;
	      keplerian_t orbit = keplerian_for_cartesian( b, center );
	      if(planets_coordinate_system==barycentric) center.mass += b.mass;
	      const double rad2deg = 180./M_PI;
	      snprintf(buf, bufsize, "%10d %lg  %6d %6d  %lg  % 9.5lg % 9.5lg % 9.5lg  % 9.5lg % 9.5lg % 9.5lg  %d", lr.msgid(), time, sys, bod, b.mass, orbit.a, orbit.e , orbit.i*rad2deg, orbit.O*rad2deg, orbit.w *rad2deg, orbit.M*rad2deg, flags);
	    }
	  if(!keplerian_output)
	    {
	      double x = b.x - center.x;
	      double y = b.y - center.y;
	      double z = b.z - center.z;
	      double vx= b.vx- center.vx;
	      double vy= b.vy- center.vy;
	      double vz= b.vz- center.vz;
	      snprintf(buf, bufsize, "%10d %lg  %6d %6d  %lg  %9.5lg %9.5lg %9.5lg  %9.5lg %9.5lg %9.5lg  %d", lr.msgid(), time, sys, bod, b.mass, x, y, z, vx, vy, vz, flags);
	    }
	  out << buf; //  << "\n";
	}
	return out;
}
    

// EVT_EJECTION
// WARNING: Logging for ejections not yet fully implemented
    std::ostream& record_output_2(std::ostream &out, gpulog::logrecord &lr, swarm::body_range_t &body_range)
{
	double T;
	int sys, body_id;
	swarm::body b;
	lr >> T >> sys >> b;

	if(!body_range.in(b.body_id)) return out;

        size_t bufsize = 1000;
        char buf[bufsize];

	if( keplerian_output) 
	  {
	    double mass = 1. + b.mass; // + center.mass;
	    keplerian_t orbit;
	    calc_keplerian_for_cartesian(orbit.a, orbit.e, orbit.i, orbit.O, orbit.w, orbit.M, b.x, b.y,b.z,b.vx,b.vy,b.vz,mass );
	    const double rad2deg = 180./M_PI;
	    snprintf(buf, bufsize, "%10d %lg  %6d %6d  %lg  % 9.5lg % 9.5lg % 9.5lg  % 9.5lg % 9.5lg % 9.5lg   ***WARNING: NOT ACCURATE***", lr.msgid(), T, sys, body_id, b.mass, orbit.a, orbit.e , orbit.i*rad2deg, orbit.O*rad2deg, orbit.w *rad2deg, orbit.M*rad2deg);
	  }
	else
	  {
	    snprintf(buf, bufsize, "%10d %lg  %6d %6d  %lg  % 9.5lf % 9.5lf % 9.5lf  % 9.5lf % 9.5lf % 9.5lf", lr.msgid(), T, sys, body_id, b.mass, b.x, b.y, b.z, b.vx, b.vy, b.vz);
	  }
	out << buf;

	return out;
}


// EVT_TRANSIT
std::ostream& record_output_15(std::ostream &out, gpulog::logrecord &lr, swarm::body_range_t &body_range)
{
	double T;
	int sys, body_id;
	double minb, vproj;
	lr >> T >> sys >> body_id >> minb >> vproj;
	if(!body_range.in(body_id)) return out;

        size_t bufsize = 1000;
        char buf[bufsize];
	snprintf(buf, bufsize, "%10d %lg  %6d %6d  %lg %lg Transit Experimental", lr.msgid(), T, sys, body_id, minb, vproj);
	out << buf;

	return out;
}


// EVT_RV_OBS
std::ostream& record_output_11(std::ostream &out, gpulog::logrecord &lr, swarm::body_range_t &body_range)
{
	double T;
	int sys, body_id;
	double vz;
	lr >> T >> sys >> body_id >> vz;
	if(!body_range.in(body_id)) return out;

        size_t bufsize = 1000;
        char buf[bufsize];
	snprintf(buf, bufsize, "%10d %lg  %6d %6d  %lg Rv Experimental", lr.msgid(), T, sys, body_id, vz);
	out << buf;

	return out;
}


// EVT_OCCULTATION
std::ostream& record_output_16(std::ostream &out, gpulog::logrecord &lr, swarm::body_range_t &body_range)
{
	double T;
	int sys, body_id;
	double minb, vproj;
	lr >> T >> sys >> body_id >> minb >> vproj;
	if(!body_range.in(body_id)) return out;

        size_t bufsize = 1000;
        char buf[bufsize];
	snprintf(buf, bufsize, "%10d %lg  %6d %6d  %lg %lg Occultation Experimental", lr.msgid(), T, sys, body_id, minb, vproj);
	out << buf;

	return out;
}


    std::ostream &output_record(std::ostream &out, gpulog::logrecord &lr, swarm::body_range_t &bod)
{
	int evtid = lr.msgid();

	// Make sure these stay synced with src/swarm/log/log.hpp
	switch(evtid){
	case 1: // standard system snapshot
	  return record_output_1(out,lr,bod);
	case 2: // data one body upon ejection
	  return record_output_2(out,lr,bod);
	case 3: // reserved for data for one pair of bodies upon close encounter/collision
	  return record_output_default(out,lr); // feature still missing
	case 11: // star v_z at observation time
	  return record_output_11(out,lr,bod);
	case 15: // near a transit of planet in front of star
	  return record_output_15(out,lr,bod);
	case 16: // near an occultation of star in front of planet
	  return record_output_16(out,lr,bod);
	default:
	  return record_output_default(out,lr);
	}
}


    //    void execute(const std::string &datafile, swarm::time_range_t T, swarm::sys_range_t sys)
      void execute(const std::string &datafile, swarm::time_range_t T, swarm::sys_range_t sys, swarm::body_range_t bod = swarm::body_range_t() )
{
	swarmdb db(datafile);
	swarmdb::result r = db.query(sys, T);
		//	swarmdb::result r = db.query(sys, bod, T);
	gpulog::logrecord lr;
	while(lr = r.next())
	{
	  output_record(std::cout, lr,bod );
		std::cout << "\n";
	}
}

  } } // end namespace swarm::query
