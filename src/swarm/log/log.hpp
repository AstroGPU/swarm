/*************************************************************************
 * Copyright (C) 2010 by Mario Juric and the Swarm-NG Development Team   *
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

#pragma once

#include "../common.hpp"

#include "gpulog/gpulog.h"
#include "gpulog/lprintf.h"

#include "../types/ensemble.hpp"

#include "types.hpp"

namespace swarm { 
/**
\title The event/data logging system for swarm

Mario Juric

Despite the horrendously complicated implementation in log.hpp, the ideas
behind it are quite simple:

* There are two "event log buffer" objects to which the user can log events,
  or the current states of systems being integrated: hlog and dlog.  They
  should be viewed as two facets of the same object: hlog the facet exposed
  towards the CPU code, while dlog is to be used by the GPU code. 
  Otherwise, their interfaces are nearly identical.  Currently, these are
  global objects, both on CPU and GPU (somewhat like cout in C++)

* There are two types of messages you can send to these log buffers. They
  differ in that they get stored to their individual "sub-buffers", each 
  being optimized for the message of a given type:

  1) Bodies - using log_body(ens, sys, bod, T) you can store the current
              (m,x,v,T) coordinates of a body bod in system sys, plus an
              integer worth of arbitrary user data. Look at 
              eventlog_base::body structure to see exactly what is stored.
              This is the basic function with which snapshots of full system
              state can be built (see docs/snapshotting.txt).

  2) Events - these are general purpose records, not longer than ~200 bytes,
              that begin with an integer eventId and then store whatever
              data you need to store.  Use log_event() set of templated
              functions to log events.  Their effective interface behaves
              like:

		\code
		void log_event(int evtId, ...)
		\endcode

              where ... stands for up to 9 variables of arbitrary types
              (think printf from C).  The machinery beneath these functions
              will byte-copy the value of all arguments into the event
              buffer.  As the name suggests, use them to record events. 

              Example usage:

		\code
                #define EVT_EJECTTION 1
                #define EVT_ENCOUTNER 2
                   ....
                dlog.log_event(EVT_ENCOUTNER, bod1, bod2);
		\endcode

  2a) Printf- an implementation of C printf() functionality, callable from
              the GPU.  Internally, these are just events of event ID
              EVT_PRINTF, but they also have a set of convenience template
              functions that makes the interface look exactly like that of
              printf().  On the CPU side, there's some code to recognize
              EVT_PRINTF events, reconstruct them and pass to C's printf()
              to produce the output string. GPU interface example:

			\code
                dlog.printf("Debug: body %d just hit body %d", bod1, bod2);
			\endcode

	      where bod1, bod2 are integers.

  From the integrator-developer's point of view, these functions allow an
  easy way to record that something has happened in the simulation, and to
  store a (sub)set of bodies related to that event of interest.  These
  events and accompanying data will be marshalled to output files for
  subsequent data analysis. For example, let's say you want to record a
  collision event between bodies bod1 and bod2. You might do it like this:

  \code
    int evtref = dlog.log_event(EVT_ENCOUTNER, sys, bod1, bod2, T);
    dlog.log_body(ens, sys, bod1, T, evtref);
    dlog.log_body(ens, sys, bod2, T, evtref);
\endcode

  The above will store a user-defined EVT_ENCOUTNER event to the event log,
  recording the system id, the ids of the two bodies involved, and the time
  of the collision.  log_event returns a unique integer (evtref). 
  Subsequent two lines stored all the information about the bodies (masses,
  position, velocities..), together with the event reference (evtref) which
  can be used in the subsequent analysis to match these two entries in the
  outputs with the reason why they're there.

* For now, the event buffers cannot be read directly but are automatically
  "piped" to a "writer" object.  A writer is an object that inherits from
  abstract class writer (see swarm.h), and is attached to hlog using hlog's
  attach_sink() function.  writer implements 'void process(ieventstream&)'
  method that when called should drain the buffers of all events (usually
  storing them to a file, or outputting some to the screen (e.g., the
  printfs)).  This method is called whenever the event buffer is flushed
  (more below).

  The argument, an ieventstream object, provides an istream-like interface to
  the eventlog.  Example: if you stored the following event (e.g., either in
  GPU code, or CPU code, doesn't matter):

	  \code
    #define EVT_MYEVENT 1
    int i = 2; double x = 3.; float2 xy = {2., 4.}
    dlog.log_event(EVT_MYEVENT, i, x, xy);
	\endcode

  when passed an ieventstream object (named 'es'), you'd read it as follows:

  \code
    switch(es.next())// advances to next event, returning its event ID
    {
    case EVT_MYEVENT:
       int i; double x; float2 xy;
       es >> i >> x >> xy;
       ... do something useful ...
       break;
    ...
    }
	\endcode

  This is only for reading events. For accessing the sub-buffer holding the
  bodies, an array-like interface is provided instead.

  Sample implementation (but still very useful) of a writer is class
  binary_writer, in swarmlog.cpp.  binary_writer just binary-dumps the
  bodies and events sub-buffer into two files, while printf-ing any
  EVT_PRINTF events to the screen.

* Both the CPU and GPU buffers have limited reserved space (which is why I
  refer to them as buffers).  Periodically, they have to be flushed;
  otherwise they'll start drooping events/bodies.

  Flushing is accomplished by calling hlog.flush() (immediately), or via
  hlog.flush_if_needed() (that will flush() only if buffers are near
  capacity).  There's no explicit call to flush the GPU buffer --
  hlog.flush* family will take care of that.

  Best practice: call flush_if_needed() after every kernel call (for GPU
  kernels), or after every step (CPU kernels).

* Finally, before they can be used the eventlog buffers must be initialized
  by calling hlog.initialize(), and a sink must be attached via
  hlog.attach_sink().  See swarm.cpp for an example. Also, immediately
  before launching a kernel from which you'll use dlog, you MUST call
  hlog.prepare_for_gpu(). Failure to do so will result in a rupture in
  space-time continuum and death of countless kittens. And you don't want to
  be a kitten-killer.

* Post-final notes: The implementation is rather dirty at this point, as I
  didn't have the time to make it prettier.  The design has passed through a
  few redesigns, and there may be vestiges of the old code here and there
  (so if something seems wrong or illogical, that may be the reason).  For
  the same reason, a few variable names may be misnomers.  There are also
  likely to be many, many, many bugs. But the code has been confirmed to
  work both on CPU and a GPU (a GTX 260).
  */
namespace log { /////////////////////////////////////////////////////////////////////////////////////////////

static const int EVT_SNAPSHOT		= 1;	//! marks a snapshot of a system. see swarm::log::system() down below
static const int EVT_EJECTION		= 2;	//! marks an ejection event
static const int EVT_ENCOUNTER		= 3;	//! marks an encounter event


template<typename L, typename T1>
GENERIC PTR_T(SCALAR(T1)) event(L &l, const int recid, const double T, const int sys, const T1 &v1)
	{
		return l.write(recid, T, sys, v1);
	}

template<typename L, typename T1, typename T2>
GENERIC PTR_T(SCALAR(T2)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2)
	{
		return l.write(recid, T, sys, v1, v2);
	}

template<typename L, typename T1, typename T2, typename T3>
GENERIC PTR_T(SCALAR(T3)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3)
	{
		return l.write(recid, T, sys, v1, v2, v3);
	}

template<typename L, typename T1, typename T2, typename T3, typename T4>
GENERIC PTR_T(SCALAR(T4)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4)
	{
		return l.write(recid, T, sys, v1, v2, v3, v4);
	}

template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5>
GENERIC PTR_T(SCALAR(T5)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5)
	{
		return l.write(recid, T, sys, v1, v2, v3, v4, v5);
	}

template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
GENERIC PTR_T(SCALAR(T6)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6)
	{
		return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6);
	}

template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
GENERIC PTR_T(SCALAR(T7)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7)
	{
		return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6, v7);
	}

template<typename L, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
GENERIC PTR_T(SCALAR(T8)) event(L &l, const int recid, const double T, const int sys, const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6, const T7 &v7, const T8 &v8)
	{
		return l.write(recid, T, sys, v1, v2, v3, v4, v5, v6, v7, v8);
	}


/**
 *	Store a snapshot of the entire system (EVT_SNAPSHOT).
 */
template<typename L>
GENERIC void system(L &l, ensemble::SystemRefConst sys)
{
	body *bodies = swarm::log::event(l, EVT_SNAPSHOT, sys.time(), sys.id() 
			, sys.state() , sys.nbod(), gpulog::array<body>(sys.nbod()));

	if(bodies != NULL) // buffer overflow hasn't happened
	{
		for(int bod=0; bod < sys.nbod(); bod++)
		{
			bodies[bod].set(bod,sys[bod]);
		}
	}
}

/**
 *	Store a snapshot of the entire system (EVT_SNAPSHOT).
 */
template<typename L>
GENERIC void system(L &l, const ensemble &ens, const int sys) { system(l,ens[sys]); }

/**
 *	Store a snapshot of the entire ensemble (convenience).
 */
template<typename L>
GENERIC void ensemble(L &l, const swarm::ensemble &ens)
{
	for(int sys = 0; sys < ens.nsys(); sys++)
		system(l, ens[sys]);
}

/**
 * Store a snapshot of all non-disabled systems within the ensemble (convenience).
 */
template<typename L>
GENERIC void ensemble_enabled(L &l, const swarm::ensemble &ens)
{
        for(int sys = 0; sys < ens.nsys(); sys++)
           {
           if(ens[sys].is_enabled())
             system(l, ens[sys]);
           }
}


} } 
