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

/*! \file device_settings.hpp
 *   \brief Declaration of funtions for GPU device setting
 *
 */

void select_cuda_device(int dev);
void print_device_information();
int blocks_per_mp( int blocksize, int shmem_per_block ) ;
bool check_cuda_limits ( int blocksize, int shmem_per_block );
int optimized_system_per_block(int chunk_size, int thread_per_system
		, int shmem_per_system, int& block_count);
void set_more_cache();
