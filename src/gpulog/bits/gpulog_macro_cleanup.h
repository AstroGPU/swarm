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

#ifndef bits_gpulog_macro_cleanup_h__
#define bits_gpulog_macro_cleanup_h__

	// Macro cleanup
	#ifdef LOCALLY_DEFINED_ALIGN
		#undef ALIGN
		#undef LOCALLY_DEFINED_ALIGN
	#endif

	#undef ISARRAY
//	#undef SCALAR
	#undef ISUNSPEC
	#undef A
	#undef M
	#undef XSTART
	#undef XSIZEOF
	#undef ASTART
	#undef SIZEOF
	#undef ADDR
//	#undef PTR_T

#endif // bits_gpulog_macro_cleanup_h__
