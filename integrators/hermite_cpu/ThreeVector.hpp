/*************************************************************************
 * Copyright (C) 2001-2010 by Eric Ford & the Swarm-NG Development Team  *
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

/*! \file ThreeVector.hpp
 *  \brief ThreeVector template for hermite_cpu
 *
 *  \todo merge various ThreeVector files into just one (or two if need separate for g++ and nvcc)
*/ 

#ifndef H_THREEVECTOR
#define H_THREEVECTOR
#include<cassert>
#include<vector>

template<typename T>
class ThreeVector
{
 public:
  typedef T FloatT;
  typedef ThreeVector<FloatT> ThreeVectorT;
 private:
#ifdef DEBUG
 public:
#endif
  FloatT mX, mY, mZ;
 public:

  ThreeVector() : mX(0.), mY(0.), mZ(0.) {};
  ThreeVector(const FloatT val) : mX(val), mY(val), mZ(val) {};
  ThreeVector(const FloatT &x, const FloatT &y, const FloatT &z=0.) : mX(x), mY(y), mZ(z) {};
  ThreeVector(const std::vector<T>& v)
    {
      assert(v.size()==3);
      mX = v[0];
      mY = v[1];
      mZ = v[2];
    };
  FloatT& X() { return mX; };
  FloatT& Y() { return mY; };
  FloatT& Z() { return mZ; };
  FloatT X() const { return mX; };
  FloatT Y() const { return mY; };
  FloatT Z() const { return mZ; };

  ThreeVector(const ThreeVectorT& v) : mX(v.X()), mY(v.Y()), mZ(v.Z()) {};
  FloatT MagnitudeSquared() const { return mX*mX+mY*mY+mZ*mZ; };
  FloatT Magnitude() const { return sqrt(MagnitudeSquared()); };
  FloatT Normalize(FloatT Mag) { mX/=Mag;  mY/=Mag;  mZ/=Mag;  return Mag; };
  FloatT Normalize() { return Normalize(Magnitude()); };
  ThreeVectorT& operator+=(const ThreeVectorT &X) 
    { mX+=X.mX; mY+=X.mY; mZ+=X.mZ; return *this; };
  ThreeVectorT& operator-=(const ThreeVectorT &X) 
    { mX-=X.mX; mY-=X.mY; mZ-=X.mZ; return *this; };
  ThreeVectorT& operator*=(const FloatT &X)
    { mX*=X; mY*=X; mZ*=X; return *this; };
  ThreeVectorT& operator/=(const FloatT &X)
    { mX/=X; mY/=X; mZ/=X; return *this; };
  ThreeVectorT& operator/=(const ThreeVectorT &X) 
    { Normalize(X.Magnitude());  return *this; };
  
  ThreeVectorT operator+(const ThreeVectorT X) const
    { ThreeVectorT Sum(*this);  Sum+=X;  return Sum; };
  ThreeVectorT operator-(const ThreeVectorT X) const
    { ThreeVectorT Diff(*this);  Diff-=X;  return Diff; };
  FloatT operator*(const ThreeVectorT X) const
    { return mX*X.mX+mY*X.mY+mZ*X.mZ; };
  ThreeVectorT operator*(const FloatT &X) const
    { ThreeVectorT Scaled(*this);  Scaled*=X;  return Scaled; };
  ThreeVectorT operator/(const FloatT &X) const
    { ThreeVectorT Scaled(*this);  Scaled/=X;  return Scaled; };
  ThreeVectorT operator-() const
    { ThreeVectorT Diff;  Diff.mX=-mX;  Diff.mY=-mY;  Diff.mZ=-mZ;  return Diff; };
  bool operator==(const ThreeVector &X) const
  { if((X()==X.X())&&(Y()==X.Y())&&(Z()==X.Z())) return true;  else return false; };
  bool Approx(const ThreeVectorT &X) const
    { if(!approx(mX,X.mX)||!approx(mY,X.mY)||!approx(mZ,X.mZ)) return 0; else return 1; };

  //  Next two lines don't work with cygwin g++ for some reason
  //  friend ThreeVectorT operator*<>(const FloatT &X, const ThreeVectorT &Y);
  //  friend ThreeVectorT cross<>(const ThreeVectorT &a, const ThreeVectorT &b);

  ThreeVectorT& RotateZ(const FloatT s, const FloatT c) throw() 
    {
      assert(approx(c*c+s*s,1.,1.e-6));
      ThreeVectorT C(*this);
      X() = + C.X()*c - C.Y()*s;
      Y() = + C.X()*s + C.Y()*c;
      return *this;
    };

  ThreeVectorT& RotateX(const FloatT s, const FloatT c)  throw() 
    {
      assert(approx(c*c+s*s,1.,1.e-6));
      ThreeVectorT C(*this);
      Y() = + C.Y()*c - C.Z()*s;
      Z() = + C.Y()*s + C.Z()*c;
      return *this;
    };

  ThreeVectorT& RotateY(const FloatT s, const FloatT c) throw() 
    {
      assert(approx(c*c+s*s,1.,1.e-6));
      ThreeVectorT C(*this);
      Z() = + C.Z()*c - C.X()*s;
      X() = + C.Z()*s + C.X()*c;
      return *this;
    };

  ThreeVectorT& RotateZ(const FloatT &Theta) throw() 
    {
      FloatT c = cos(Theta), s = sin(Theta);
      return RotateZ(s,c);
      /*
      ThreeVectorT C(*this);
      X() = + C.X()*c - C.Y()*s;
      Y() = + C.X()*s + C.Y()*c;
      return *this;
      */
    };

  ThreeVectorT& RotateX(const FloatT &Theta)  throw() 
    {
      FloatT c = cos(Theta), s = sin(Theta);
      return RotateX(s,c);
      /*
      ThreeVectorT C(*this);
      Y() = + C.Y()*c - C.Z()*s;
      Z() = + C.Y()*s + C.Z()*c;
      return *this;
      */
    };

  ThreeVectorT& RotateY(const FloatT &Theta) throw() 
    {
      FloatT c = cos(Theta), s = sin(Theta);
      return RotateY(s,c);
      /*
      ThreeVectorT C(*this);
      Z() = + C.Z()*c - C.X()*s;
      X() = + C.Z()*s + C.X()*c;
      return *this;
      */
    };

  ThreeVectorT& Rotate(const FloatT &Theta, const ThreeVectorT &Axis ) throw() 
    {
      FloatT c = cos(Theta), s = sin(Theta);
      ThreeVectorT V = Axis;
      if(V.Magnitude()==0.)  return *this;
      V.Normalize();
      ThreeVectorT C(*this);
      X() = C.X()*(V.X()*V.X()*(1.-c)+     c) + C.Y()*(V.X()*V.Y()*(1.-c)-V.Z()*s) + C.Z()*(V.X()*V.Z()*(1.-c)+V.Y()*s);
      Y() = C.X()*(V.X()*V.Y()*(1.-c)+V.Z()*s) + C.Y()*(V.Y()*V.Y()*(1.-c)+     c) + C.Z()*(V.Y()*V.Z()*(1.-c)-V.X()*s);
      Z() = C.X()*(V.X()*V.Z()*(1.-c)-V.Y()*s) + C.Y()*(V.Y()*V.Z()*(1.-c)+V.X()*s) + C.Z()*(V.Z()*V.Z()*(1.-c)+     c);
      return *this;
    };

  FloatT ProjectOntoAxis(const ThreeVectorT& axis) const
    {
      FloatT val = *this * axis;
      val /= axis.Magnitude();
      return val;
    };

  ThreeVectorT ProjectPerpToAxis(const ThreeVectorT& Axis)  const
    {
      FloatT projfac = ProjectOntoAxis(Axis)/Axis.Magnitude();
      ThreeVectorT proj = *this;
      proj -= (Axis*projfac);
      return proj;
    };
};


template<typename T>
inline ThreeVector<T> operator*(const T &X, const ThreeVector<T> &Y)
{ ThreeVector<T> Scaled(Y);  Scaled*=X;  return Scaled; };


template<typename T>
inline T dot(const ThreeVector<T> &a, const ThreeVector<T> &b)
{
T c = a.X()*b.X()+a.Y()*b.Y()+a.Z()*b.Z();
return c;
}

template<typename T>
inline ThreeVector<T> cross(const ThreeVector<T> &a, const ThreeVector<T> &b)
{
ThreeVector<T> c;
c.X() = a.Y()*b.Z() - a.Z()*b.Y();
c.Y() = a.Z()*b.X() - a.X()*b.Z();
c.Z() = a.X()*b.Y() - a.Y()*b.X();
return c;
}

#include <iostream>
template<typename T>
std::ostream& operator<<(std::ostream& ios, const ThreeVector<T>& v)
{
  ios << v.X() << " " << v.Y() << " " << v.Z();
  return ios;
}

#endif
