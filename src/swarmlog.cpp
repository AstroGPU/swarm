#include "swarmlog.h"
#include <iostream>
#include <astro/binarystream.h>
#include <astro/memorymap.h>
#include <fstream>
#include <sstream>
#include <memory>
#include <limits>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

extern "C" void debug_hook()
{
	// a hook into which to place a breakpoint when the debugger
	// won't stop in nvcc-compiled code.
	std::cerr << "";
}

// global pointer to log output writer
namespace swarm
{
	namespace log
	{
		std::auto_ptr<writer> log_writer;
	}
}

void swarm::log::init(const std::string &writer_cfg, int host_buffer_size, int device_buffer_size)
{
	log_writer.reset(writer::create(writer_cfg));

	// log memory allocation
	hlog.alloc(host_buffer_size);
	gpulog::alloc_device_log("dlog", device_buffer_size);
}

void swarm::log::shutdown()
{
	swarm::log::init("null");
}

void swarm::log::flush(int flags)
{
	// TODO: Implement flushing of writer as well
	assert(flags & memory);

	// TODO: Implement flush-only-if-near-capacity (observe the if_full flag)
	if(!log_writer.get())
	{
		std::cerr << "No output writer attached!\n";
		assert(0);
	}

	// flush the CPU and GPU buffers
	replay_printf(std::cerr, hlog);
	log_writer->process(hlog.internal_buffer(), hlog.size());

	copy(hlog, "dlog", gpulog::LOG_DEVCLEAR);
	replay_printf(std::cerr, hlog);
	log_writer->process(hlog.internal_buffer(), hlog.size());

	hlog.clear();
}
