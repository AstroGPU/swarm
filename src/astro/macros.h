/***************************************************************************
 *   Copyright (C) 2005 by Mario Juric   *
 *   mjuric@astro.Princeton.EDU   *
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

#ifndef _astro_macros_h
#define _astro_macros_h

#define FOREACH2j(C, i_, x) for(C i_ = (x).begin(); i_ != (x).end(); ++i_)
#define FOREACH2(C, x) FOREACH2j(C, i, x)

#define REVEACH2j(C, i_, x) for(C i_ = (x).rbegin(); i_ != (x).rend(); ++i_)
#define REVEACH2(C, x) REVEACH2j(C, i, x)

#define FOREACHj(i_, x) for(typeof((x).begin()) i_ = (x).begin(); i_ != (x).end(); ++i_)
#define FOREACH(x) FOREACHj(i, x)

#define REVEACHj(i_, x) for(typeof((x).rbegin()) i_ = (x).rbegin(); i_ != (x).rend(); ++i_)
#define REVEACH(x) REVEACHj(i, x)

#define FORj(i, i0, i1) for(int i = i0; i != i1; ++i)
#define FOR(i0, i1) FORj(i, i0, i1)

#define REVj(i, i0, i1) for(int i = i0; i != i1; --i)
#define REV(i0, i1) REVj(i, i0, i1)

#define OSTREAM(T) std::ostream &operator <<(std::ostream &out, T)
#define ISTREAM(T) std::istream &operator >>(std::istream &in, T)

#endif // _astro_macros_h
