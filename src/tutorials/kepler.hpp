/*************************************************************************
 * Copyright (C) 2009-2010 by Eric Ford & the Swarm-NG Development Team  *
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


/*! \file kepler.hpp
 *  \brief converts between Cartesian & Keplerian coordinates
 *
 *  Based on code from John Chambers' Mercury code
 *
 *  Cleaned up by Saleh on October 11, 2011. First I split this
 *  up into Header/Source files and used it in montecarlo.cpp
*/
#pragma once

double improve_mean_to_eccentric_annomaly_guess(const double e, const double M, const double x);
double mean_to_eccentric_annomaly(const double e,  double M);
void calc_cartesian_for_ellipse(double& x,double& y, double & z, double &vx, double &vy, double &vz, const double a, const double e, const double i, const double O, const double w, const double M, const double GM);
void calc_keplerian_for_cartesian( double& a,  double& e,  double& i,  double& O,  double& w,  double& M, const double x,const double y, const double z, const double vx, const double vy, const double vz, const double GM);




