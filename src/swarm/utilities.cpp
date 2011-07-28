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

/*! \file utilities.cpp
 *  \brief provides several utility  functions for public interface for swarm libaray
 *
*/

#include "utilities.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include "io.hpp"
//
// Utilities
//


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

        return true;
}


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


}
