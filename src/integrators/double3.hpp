
#include "vector_types.h"
#include "vector_functions.h"

GENERIC double3 operator+(const double3 &a,const double3 &b) {
  return make_double3(a.x+b.x, a.y+b.y, a.z+b.z);
}
GENERIC double3 operator-(const double3 &a,const double3 &b) {
  return make_double3(a.x-b.x, a.y-b.y, a.z-b.z);
}
GENERIC double3 operator*(const double3 &a, const double& b) {
  return make_double3(a.x*b, a.y*b, a.z*b);
}
GENERIC double3 operator*(const double& b, const double3& a) {
  return a * b;
}
GENERIC double3 operator/(const double3 &a, const double& b){
  return a *(1.0/b);
}

GENERIC double3& operator+=(double3& a, const double3& b){
  a = a + b;
}
GENERIC double3& operator-=(double3& a, const double3& b){
  a = a - b;
}
 

GENERIC double sqnorm(double3 a) { return a.x*a.x+a.y*a.y+a.z*a.z; }
GENERIC double inner_product(double3 a,double3 b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
