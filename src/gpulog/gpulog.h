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
