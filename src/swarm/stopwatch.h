/*************************************************************************
 * Copyright (C) 2010 by Mario Juric  and the Swarm-NG Development Team  *
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

/*! \file stopwatch.h
 * \brief Defines stopwatch class for benchmarking cpu & gpu performance
 * 
 * class stopwath is based on NVIDIA's LinuxStopWatch class
 */

#ifndef stopwatch_h__
#define stopwatch_h__

#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdint.h>

/// Class for benchmarking cpu & gpu performance.  Based on NVIDIA's LinuxStopWatch class
class stopwatch
{
protected:
	struct timeval  start_time;	// Start of measurement
	float  diff_time;		// Time difference between the last start and stop (in ms)
	float  total_time;		// TOTAL time difference between starts and stops (in ms)
	bool running;			// flag if the stop watch is running
	int clock_sessions;		// Number of times clock has been started and stopped (for averaging)

public:
	stopwatch() :
		start_time(),
		diff_time(0.0),
		total_time(0.0),
		running(false),
		clock_sessions(0)
	{ }

	/// Start time measurement
	void start()
	{
		gettimeofday( &start_time, 0);
		running = true;
	}

	/// Stop time measurement and increment add to the current diff_time summation
	/// variable. Also increment the number of times this clock has been run.
	void stop()
	{
		diff_time = getDiffTime();
		total_time += diff_time;
		running = false;
		clock_sessions++;
	}

	// Reset the timer to 0. Does not change the timer running state but does
	// recapture this point in time as the current start time if it is running.
	void reset()
	{
		diff_time = 0;
		total_time = 0;
		clock_sessions = 0;
		if( running )
		{
			gettimeofday( &start_time, 0);
		}
	}

	// Time in sec. after start. If the stop watch is still running (i.e. there
	// was no call to stop()) then the elapsed time is returned added to the
	// current diff_time sum, otherwise the current summed time difference alone
	// is returned.
	float getTime() const
	{
		// Return the TOTAL time to date
		float retval = total_time;
		if(running)
		{
			retval += getDiffTime();
		}
		return 0.001 * retval;
	}

	// Add dt seconds of time to the counter. Does not increment the number of
	// clock sessions
	void addTime(const float dt)
	{
		total_time += 1000*dt;
	}

	// Time in msec. for a single run based on the total number of COMPLETED runs
	// and the total time.
	float getAverageTime() const
	{
		return 0.001 * total_time/clock_sessions;
	}

	int nSessions() const
	{
		return clock_sessions;
	}

private:

	// helpers functions

	float getDiffTime() const
	{
		struct timeval t_time;
		gettimeofday( &t_time, 0);

		// time difference in milli-seconds
		return  (float) (1000.0 * ( t_time.tv_sec - start_time.tv_sec)
				+ (0.001 * (t_time.tv_usec - start_time.tv_usec)) );
	}
};

template< class T >
double watch_time(const T& a){
	stopwatch s;
	s.start();
	a();
	SYNC;
	s.stop();
	return s.getTime()*1000.;
}

#endif // stopwatch_h__
