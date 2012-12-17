/***************************************************************************
 *   Copyright (C) 2005 by Mario Juric                                     *
 *   mjuric@astro.Princeton.EDU                                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
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

/*! \file fakemmap.h
 *   \brief Defines macros and functions for faking memory mapping for portability.
 * Fake memory mapping for portability with windows system
 *
 *  *EXPERIMENTAL*: This module is not thoroughly tested.
 *  \ingroup experimental
 * 
 *
 *
 */


#include <stdlib.h>

#ifdef _WIN32
#include <wchar.h>
#include <io.h>
#include <fcntl.h>
inline int getpagesize() { return 4096; }
#else
#include <unistd.h>
#endif


void *fakemmap(void *addr, size_t length, int prot, int flags, int fd, off_t offset);
int fakemunmap(void *addr, size_t length);
int fakemsync(void *addr, size_t length, int flags);


#define PROT_READ       0x1             /* Page can be read.  */
#define PROT_WRITE      0x2             /* Page can be written.  */
#define PROT_EXEC       0x4             /* Page can be executed.  */
#define PROT_NONE       0x0             /* Page can not be accessed.  */
#define MAP_FAILED      ((void *) -1)
#define MAP_SHARED      0x01            /* Share changes.  */
#define MAP_PRIVATE     0x02            /* Changes are private.  */
/* Flags to `msync'.  */
#define MS_ASYNC        1               /* Sync memory asynchronously.  */
#define MS_SYNC         4               /* Synchronous memory sync.  */
#define MS_INVALIDATE   2               /* Invalidate the caches.  */

