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

  __device__ __host__ ThreeVector() : mX(0.), mY(0.), mZ(0.) {};
  __device__ __host__ ThreeVector(const FloatT val) : mX(val), mY(val), mZ(val) {};
  __device__ __host__ ThreeVector(const FloatT &x, const FloatT &y, const FloatT &z=0.) : mX(x), mY(y), mZ(z) {};
  ThreeVector(const std::vector<T>& v)
    {
//      assert(v.size()==3);
      mX = v[0];
      mY = v[1];
      mZ = v[2];
    };
  __device__ __host__ FloatT& X() { return mX; };
  __device__ __host__ FloatT& Y() { return mY; };
  __device__ __host__ FloatT& Z() { return mZ; };
  __device__ __host__ FloatT X() const { return mX; };
  __device__ __host__ FloatT Y() const { return mY; };
  __device__ __host__ FloatT Z() const { return mZ; };

  __device__ __host__ ThreeVector(const ThreeVectorT& v) : mX(v.X()), mY(v.Y()), mZ(v.Z()) {};
  __device__ __host__ FloatT MagnitudeSquared() const { return mX*mX+mY*mY+mZ*mZ; };
  __device__ __host__ FloatT Magnitude() const { return sqrt(MagnitudeSquared()); };
  __device__ __host__ FloatT Normalize(FloatT Mag) { mX/=Mag;  mY/=Mag;  mZ/=Mag;  return Mag; };
  __device__ __host__ FloatT Normalize() { return Normalize(Magnitude()); };
  __device__ __host__ ThreeVectorT& operator+=(const ThreeVectorT &X) 
    { mX+=X.mX; mY+=X.mY; mZ+=X.mZ; return *this; };
  __device__ __host__ ThreeVectorT& operator-=(const ThreeVectorT &X) 
    { mX-=X.mX; mY-=X.mY; mZ-=X.mZ; return *this; };
  __device__ __host__ ThreeVectorT& operator*=(const FloatT &X)
    { mX*=X; mY*=X; mZ*=X; return *this; };
  __device__ __host__ ThreeVectorT& operator/=(const FloatT &X)
    { mX/=X; mY/=X; mZ/=X; return *this; };
  __device__ __host__ ThreeVectorT& operator/=(const ThreeVectorT &X) 
    { Normalize(X.Magnitude());  return *this; };
  
  __device__ __host__ ThreeVectorT operator+(const ThreeVectorT X) const
    { ThreeVectorT Sum(*this);  Sum+=X;  return Sum; };
  __device__ __host__ ThreeVectorT operator-(const ThreeVectorT X) const
    { ThreeVectorT Diff(*this);  Diff-=X;  return Diff; };
  __device__ __host__ FloatT operator*(const ThreeVectorT X) const
    { return mX*X.mX+mY*X.mY+mZ*X.mZ; };
  __device__ __host__ ThreeVectorT operator*(const FloatT &X) const
    { ThreeVectorT Scaled(*this);  Scaled*=X;  return Scaled; };
  __device__ __host__ ThreeVectorT operator/(const FloatT &X) const
    { ThreeVectorT Scaled(*this);  Scaled/=X;  return Scaled; };
  __device__ __host__ ThreeVectorT operator-() const
    { ThreeVectorT Diff;  Diff.mX=-mX;  Diff.mY=-mY;  Diff.mZ=-mZ;  return Diff; };
  __device__ __host__ bool operator==(const ThreeVector &X) const
  { if((X()==X.X())&&(Y()==X.Y())&&(Z()==X.Z())) return true;  else return false; };
  __device__ __host__ bool Approx(const ThreeVectorT &X) const
    { if(!approx(mX,X.mX)||!approx(mY,X.mY)||!approx(mZ,X.mZ)) return 0; else return 1; };

  //  Next two lines don't work with cygwin g++ for some reason
  //  friend ThreeVectorT operator*<>(const FloatT &X, const ThreeVectorT &Y);
  //  friend ThreeVectorT cross<>(const ThreeVectorT &a, const ThreeVectorT &b);

  __device__ __host__ ThreeVectorT& RotateZ(const FloatT s, const FloatT c) throw() 
    {
//      assert(approx(c*c+s*s,1.,1.e-6));
      ThreeVectorT C(*this);
      X() = + C.X()*c - C.Y()*s;
      Y() = + C.X()*s + C.Y()*c;
      return *this;
    };

  __device__ __host__ ThreeVectorT& RotateX(const FloatT s, const FloatT c)  throw() 
    {
//      assert(approx(c*c+s*s,1.,1.e-6));
      ThreeVectorT C(*this);
      Y() = + C.Y()*c - C.Z()*s;
      Z() = + C.Y()*s + C.Z()*c;
      return *this;
    };

  __device__ __host__ ThreeVectorT& RotateY(const FloatT s, const FloatT c) throw() 
    {
//      assert(approx(c*c+s*s,1.,1.e-6));
      ThreeVectorT C(*this);
      Z() = + C.Z()*c - C.X()*s;
      X() = + C.Z()*s + C.X()*c;
      return *this;
    };

  __device__ __host__ ThreeVectorT& RotateZ(const FloatT &Theta) throw() 
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

  __device__ __host__ ThreeVectorT& RotateX(const FloatT &Theta)  throw() 
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

  __device__ __host__ ThreeVectorT& RotateY(const FloatT &Theta) throw() 
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

  __device__ __host__ ThreeVectorT& Rotate(const FloatT &Theta, const ThreeVectorT &Axis ) throw() 
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

  __device__ __host__ FloatT ProjectOntoAxis(const ThreeVectorT& axis) const
    {
      FloatT val = *this * axis;
      val /= axis.Magnitude();
      return val;
    };

  __device__ __host__ ThreeVectorT ProjectPerpToAxis(const ThreeVectorT& Axis)  const
    {
      FloatT projfac = ProjectOntoAxis(Axis)/Axis.Magnitude();
      ThreeVectorT proj = *this;
      proj -= (Axis*projfac);
      return proj;
    };
};


template<typename T>
__device__ __host__ inline ThreeVector<T> operator*(const T &X, const ThreeVector<T> &Y)
{ ThreeVector<T> Scaled(Y);  Scaled*=X;  return Scaled; };


template<typename T>
__device__ __host__ inline T dot(const ThreeVector<T> &a, const ThreeVector<T> &b)
{
T c = a.X()*b.X()+a.Y()*b.Y()+a.Z()*b.Z();
return c;
}

template<typename T>
__device__ __host__ inline ThreeVector<T> cross(const ThreeVector<T> &a, const ThreeVector<T> &b)
{
ThreeVector<T> c;
c.X() = a.Y()*b.Z() - a.Z()*b.Y();
c.Y() = a.Z()*b.X() - a.X()*b.Z();
c.Z() = a.X()*b.Y() - a.Y()*b.X();
return c;
}

#if !__CUDACC__	// workaround for CUDA 2.2 devemu bug
#include <iostream>
template<typename T>
std::ostream& operator<<(std::ostream& ios, const ThreeVector<T>& v)
{
  ios << v.X() << " " << v.Y() << " " << v.Z();
  return ios;
}
#endif

#endif
