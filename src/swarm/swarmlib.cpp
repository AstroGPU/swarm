/*************************************************************************
 * Copyright (C) 2008-2010 by Mario Juric & Swarm-NG Development Team    *
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

/*! \file swarmlib.cpp
 *  \brief provides several functions for public interface for swarm libaray
 *
 *  If you want to use the swarm library and have questions about what the 
 *  functions do, this is a great plae to start looking.
*/

#include <cuda_runtime_api.h>
#include "swarm.h"
#include <vector>
#include <algorithm> // for swap
#include <memory>
#include <iostream>
#include <dlfcn.h>
#include <sstream>
#include <fstream>
#include <valarray>
#include "io.hpp"
#include "cux/cux.h"
//
// Utilities
//

static bool swarm_initialized = false;

namespace swarm {

/*!
   \brief Configur grid

   Given a total number of threads, their memory requirements, and the
   number of threadsPerBlock, compute the optimal allowable grid dimensions.
   Returns false if the requested number of threads are impossible to fit to
   shared memory.

  @param[out] gridDim
  @param[in] threadsPerBlock
  @param[in] nthreads
  @param[in] dynShmemPerThread
  @param[in] staticShmemPerBlcok
  @return boolean
 */
bool configure_grid(dim3 &gridDim, int threadsPerBlock, int nthreads, int dynShmemPerThread, int staticShmemPerBlock)
{
        const int shmemPerMP =  16384;

        int dyn_shared_mem_required = dynShmemPerThread*threadsPerBlock;
        int shared_mem_required = staticShmemPerBlock + dyn_shared_mem_required;
        if(shared_mem_required > shmemPerMP) { return false; }

        // calculate the total number of threads
        int nthreadsEx = nthreads;
        int over = nthreads % threadsPerBlock;
        if(over) { nthreadsEx += threadsPerBlock - over; } // round up to multiple of threadsPerBlock

        // calculate the number of blocks
        int nblocks = nthreadsEx / threadsPerBlock;
        if(nthreadsEx % threadsPerBlock) { nblocks++; }

        // calculate block dimensions so that there are as close to nblocks blocks as possible
        find_best_factorization(gridDim.x, gridDim.y, nblocks);
        gridDim.z = 1;

/* Commenting out since never gets executed right now -- CW 8/9/10
 * #if 0
 *      std::cerr << "+ Grid configuration =========================\n";
 *      std::cerr << "      Threads requested = " << nthreads << " with shmem/thread = " << dynShmemPerThread << " and shmem/blk = " << staticShmemPerBlock << "\n";
 *      std::cerr << "      Grid configured as (" << gridDim.x << ", " << gridDim.y << ", " << gridDim.z <<") array of blocks with " << threadsPerBlock << " threads per block.\n";
 *      std::cerr << "      Total threads to execute = " << nthreadsEx << "\n";
 *      std::cerr << "- Grid configuration =========================\n";
 * #else
 *      //std::cerr << "Kernel exec. config: (" << gridDim.x << ", " << gridDim.y << ", " << gridDim.z <<") x " << threadsPerBlock << " thr/blk (" << nthreadsEx << " thr total; " << nthreads << " thr needed)\n";
 * #endif
 */
        return true;
}

//----------------------------------------------------------------------------

void debugger_stop()
{
        std::cerr << "Block for debugger here!\n";
}

//----------------------------------------------------------------------------

/*!
   \brief Find best factorization

   Find the dimensions (bx,by) of a 2D grid of blocks that has as close to nblocks blocks as possible
  @param[out] bx
  @param[out] by
  @param[in] nblocks
*/
void find_best_factorization(unsigned int &bx, unsigned int &by, int nblocks)
{
        bx = -1;
        int best_r = 100000;
        for(int bytmp = 1; bytmp != 65536; bytmp++)
        {
                int r  = nblocks % bytmp;
                if(r < best_r && nblocks / bytmp < 65535)
                {
                        by = bytmp;
                        bx = nblocks / bytmp;
                        best_r = r;

                        if(r == 0) { break; }
                        bx++;
                }
        }
        if(bx == -1) { std::cerr << "Unfactorizable?!\n"; exit(-1); }
}

//----------------------------------------------------------------------------

/*!
   \brief Initialize the swarm library.

   This function must be called before any other.
   @param[in] cfg configuration class
*/
void init(const config &cfg)
{
        if(swarm_initialized) { return; }

        // Initialize appropriate GPU
        cux_init();

        swarm_initialized = true;
}

//----------------------------------------------------------------------------

/*!
   \brief  load a configuration file

   @param[out] cfg configuration class
   @param[in] fn file name sting
 */
void load_config(config &cfg, const std::string &fn)
{
	std::ifstream in(fn.c_str());
	if(!in) ERROR("Cannot open configuration file '" + fn + "'.");

	std::string line;
	int iline = 0;
	while(std::getline(in, line))
	{
		iline++;
		line = trim(line);
		if(line.empty()) { continue; }
		if(line[0] == '#') { continue; }

		size_t eqpos = line.find('=');
		if(eqpos == std::string::npos) ERROR("Error on line " + line + ": '=' sign expected.");

		std::string key = trim(line.substr(0, eqpos)), val = trim(line.substr(eqpos+1));

		cfg[key] = val;
	}
}

//----------------------------------------------------------------------------

}
