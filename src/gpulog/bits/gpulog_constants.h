#ifndef bits_gpulog_constants_h__
#define bits_gpulog_constants_h__

namespace gpulog
{
	//
	// System message IDs
	//
	static const int MSG_INVALID = -1;
	static const int MSG_PRINTF = -2;

	//
	// flags for gpulog::copy() and related functions
	//
	static const int LOG_DEVCLEAR = 0x01;
}

#endif // bits_gpulog_constants_h__
