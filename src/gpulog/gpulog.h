/***************************************************************************
 *   Copyright (C) 2010 by Mario Juric   *
 *   mjuric@cfa.harvard.EDU       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
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

#ifndef gpulog_h__
#define gpulog_h__

#include <cstdio>
#include <cassert>

#ifndef ARGINFO
	#define ARGINFO 0
#endif

#if ARGINFO
	#define IFARGINFO(x) x
#else
	#define IFARGINFO(x)
#endif

#include "bits/gpulog_debug.h"
#include "bits/gpulog_align.h"
#include "bits/gpulog_types.h"
#include "bits/gpulog_ttraits.h"
#include "bits/gpulog_constants.h"
#include "bits/gpulog_msg_layout.h"
#include "bits/gpulog_log.h"
#include "bits/gpulog_logrecord.h"
#include "bits/gpulog_ilogstream.h"
#include "bits/gpulog_macro_cleanup.h"

namespace gpulog
{
	// Import the externally visible classes and constant
	// into gpulog namespace
	using internal::device_log;
	using internal::host_log;
	using internal::array;
	using internal::header;
	using internal::logrecord;
	using internal::ilogstream;

	using internal::alloc_device_log;
	using internal::free_device_log;
}

#endif // gpulog_h__
