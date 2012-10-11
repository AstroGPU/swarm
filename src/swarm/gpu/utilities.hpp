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

/*! \file utilities.hpp
 *   \brief Defines a function for counting the number of active ensemble 
 *          systems on GPU and utility functions for finding best factorization
 *          and compute the optimal allowable grid dimensions based on 
 *          the system configuration. 
 */

#pragma once
#include "../types/ensemble.hpp"

namespace swarm {

/*!
 *  \brief Count the number of active systems on GPU without downloading the ensemble.
 *  It is a very small kernel that iterates through all systems and counts
 *  The number of active ones. One can optimize the performance of this using
 *  Prefix-sum algorithm.
 *
 *  But this is such a small kernel and it is not very frequently used.
 */
int number_of_active_systems(deviceEnsemble ens);

void reactivate_systems(deviceEnsemble ens);

/*!
   \brief Find best factorization

   Find the dimensions (bx,by) of a 2D grid of blocks that has as close to nblocks blocks as possible
  @param[out] bx
  @param[out] by
  @param[in] nblocks
*/
void find_best_factorization(unsigned int &bx, unsigned int &by, int nblocks);

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
bool configure_grid(dim3 &gridDim, int threadsPerBlock, int nthreads, int dynShmemPerThread, int staticShmemPerBlock);
}
