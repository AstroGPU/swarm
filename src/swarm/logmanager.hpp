/*************************************************************************
 * Copyright (C) 2011 by Saleh Dindar and the Swarm-NG Development Team  *
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
#include "gpulog/gpulog.h"
#include "gpulog/lprintf.h"
#include "log.hpp"
#include "writer.h"

#include <boost/shared_ptr.hpp>

namespace swarm { namespace log {

	class manager {
		gpulog::host_log hlog;
		gpulog::device_log* pdlog;
		std::auto_ptr<writer> log_writer;
		public:

		enum { memory = 0x01, if_full = 0x02 };

		void init(const config&, int host_buffer_size = 50*1024*1024, int device_buffer_size = 50*1024*1024);
		void flush(int flags = memory);
		void shutdown();

		gpulog::device_log* get_gpulog() { return pdlog; }
		gpulog::host_log* get_hostlog() { return &hlog; }

	};

}}
