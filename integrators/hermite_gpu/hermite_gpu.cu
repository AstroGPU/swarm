#include "swarm.h"
#include "hermite_gpu.h"

namespace gpu_hermite_aux
{
	//
	// Wrap all aux. functions in a separate namespace, to avoid
	// collisions with equally named functions from other integrators.
	//

	__device__ float3 operator*(const float3 &a, const float &b)
	{
		return make_float3(a.x*b, a.y*b, a.z*b);
	}

	__device__ float3 operator+(const float3 &a, const float3 &b)
	{
		return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
	}

}

//#if PRECISION == 3
// Double precision
#define RSQRT(x) rsqrtf(x)
#define SQRT(x)   sqrtf(x)
//typedef double real;
//#else
//// Single or Mixed precision
//#define RSQRT(x) rsqrtf(x)
//#define SQRT(x)   sqrtf(x)
//typedef float real;
//#endif

inline __device__ void copyArray(real *target, real *source, int numArray)
{
	if(numArray==9){
		target[0]=source[0];
		target[1]=source[1];
		target[2]=source[2];
		target[3]=source[3];
		target[4]=source[4];
		target[5]=source[5];
		target[6]=source[6];
		target[7]=source[7];
		target[8]=source[8];
	}
	else {
		for(int i=0; i<numArray; i++)
			target[i]=source[i];
	}
}

inline __device__ void doubleTofloat(float *floatA, real *doubleA, int numArray)
{
	if(numArray==9){
		floatA[0]=__double2float_rn(doubleA[0]);
		floatA[1]=__double2float_rn(doubleA[1]);
		floatA[2]=__double2float_rn(doubleA[2]);
		floatA[3]=__double2float_rn(doubleA[3]);
		floatA[4]=__double2float_rn(doubleA[4]);
		floatA[5]=__double2float_rn(doubleA[5]);
		floatA[6]=__double2float_rn(doubleA[6]);
		floatA[7]=__double2float_rn(doubleA[7]);
		floatA[8]=__double2float_rn(doubleA[8]);
	}
	else {
		for(int i=0; i<numArray; i++)
			floatA[i]=__double2float_rn(doubleA[i]);
	}
}

inline __device__ void floatTodouble(double *doubleA, float *floatA, int numArray)
{
	if(numArray==9){
		doubleA[0]=(double)floatA[0];
		doubleA[1]=(double)floatA[1];
		doubleA[2]=(double)floatA[2];
		doubleA[3]=(double)floatA[3];
		doubleA[4]=(double)floatA[4];
		doubleA[5]=(double)floatA[5];
		doubleA[6]=(double)floatA[6];
		doubleA[7]=(double)floatA[7];
		doubleA[8]=(double)floatA[8];
	}
	else {
		for(int i=0; i<numArray; i++)
			doubleA[i]=(double)floatA[i];
	}
}

inline __device__ void predict(real *mPos, real *mVel, real *mAcc, real *mJerk, const real dtby2, const real dtby3, double h, int numArray)
{
	if(numArray==9){
		mPos[0] += h* (mVel[0]+ dtby2*(mAcc[0]+dtby3*mJerk[0]));
		mPos[1] += h* (mVel[1]+ dtby2*(mAcc[1]+dtby3*mJerk[1]));
		mPos[2] += h* (mVel[2]+ dtby2*(mAcc[2]+dtby3*mJerk[2]));
		mPos[3] += h* (mVel[3]+ dtby2*(mAcc[3]+dtby3*mJerk[3]));
		mPos[4] += h* (mVel[4]+ dtby2*(mAcc[4]+dtby3*mJerk[4]));
		mPos[5] += h* (mVel[5]+ dtby2*(mAcc[5]+dtby3*mJerk[5]));
		mPos[6] += h* (mVel[6]+ dtby2*(mAcc[6]+dtby3*mJerk[6]));
		mPos[7] += h* (mVel[7]+ dtby2*(mAcc[7]+dtby3*mJerk[7]));
		mPos[8] += h* (mVel[8]+ dtby2*(mAcc[8]+dtby3*mJerk[8]));
		mVel[0] += h* (mAcc[0]+ dtby2*mJerk[0]);
		mVel[1] += h* (mAcc[1]+ dtby2*mJerk[1]);
		mVel[2] += h* (mAcc[2]+ dtby2*mJerk[2]);
		mVel[3] += h* (mAcc[3]+ dtby2*mJerk[3]);
		mVel[4] += h* (mAcc[4]+ dtby2*mJerk[4]);
		mVel[5] += h* (mAcc[5]+ dtby2*mJerk[5]);
		mVel[6] += h* (mAcc[6]+ dtby2*mJerk[6]);
		mVel[7] += h* (mAcc[7]+ dtby2*mJerk[7]);
		mVel[8] += h* (mAcc[8]+ dtby2*mJerk[8]);
	}
	else {
		for(int i=0; i<numArray; i++) {
			mPos[i] += h* (mVel[i]+ dtby2*(mAcc[i]+dtby3*mJerk[i]));
			mVel[i] += h* (mAcc[i]+ dtby2*mJerk[i]);
		}
	}
}

inline __device__ void correct(real *mPos, real *mVel, real *mAcc, real *mJerk, 
		real *mPosOld, real *mVelOld, real *mAccOld, real *mJerkOld, 
		const real dtby2, const real dtby6, const real dtby7, const real dt7by30, int numArray)
{
	if(numArray==9){
		mVel[0] = mVelOld[0] + dtby2*((mAccOld[0]+mAcc[0]) + dtby6*  (mJerkOld[0]-mJerk[0]));
		mVel[1] = mVelOld[1] + dtby2*((mAccOld[1]+mAcc[1]) + dtby6*  (mJerkOld[1]-mJerk[1]));
		mVel[2] = mVelOld[2] + dtby2*((mAccOld[2]+mAcc[2]) + dtby6*  (mJerkOld[2]-mJerk[2]));
		mVel[3] = mVelOld[3] + dtby2*((mAccOld[3]+mAcc[3]) + dtby6*  (mJerkOld[3]-mJerk[3]));
		mVel[4] = mVelOld[4] + dtby2*((mAccOld[4]+mAcc[4]) + dtby6*  (mJerkOld[4]-mJerk[4]));
		mVel[5] = mVelOld[5] + dtby2*((mAccOld[5]+mAcc[5]) + dtby6*  (mJerkOld[5]-mJerk[5]));
		mVel[6] = mVelOld[6] + dtby2*((mAccOld[6]+mAcc[6]) + dtby6*  (mJerkOld[6]-mJerk[6]));
		mVel[7] = mVelOld[7] + dtby2*((mAccOld[7]+mAcc[7]) + dtby6*  (mJerkOld[7]-mJerk[7]));
		mVel[8] = mVelOld[8] + dtby2*((mAccOld[8]+mAcc[8]) + dtby6*  (mJerkOld[8]-mJerk[8]));
		mPos[0] = mPosOld[0] + dtby2*((mVelOld[0]+mVel[0]) + dt7by30*((mAccOld[0]- mAcc[0]) + dtby7*(mJerkOld[0]+mJerk[0])));
		mPos[1] = mPosOld[1] + dtby2*((mVelOld[1]+mVel[1]) + dt7by30*((mAccOld[1]- mAcc[1]) + dtby7*(mJerkOld[1]+mJerk[1])));
		mPos[2] = mPosOld[2] + dtby2*((mVelOld[2]+mVel[2]) + dt7by30*((mAccOld[2]- mAcc[2]) + dtby7*(mJerkOld[2]+mJerk[2])));
		mPos[3] = mPosOld[3] + dtby2*((mVelOld[3]+mVel[3]) + dt7by30*((mAccOld[3]- mAcc[3]) + dtby7*(mJerkOld[3]+mJerk[3])));
		mPos[4] = mPosOld[4] + dtby2*((mVelOld[4]+mVel[4]) + dt7by30*((mAccOld[4]- mAcc[4]) + dtby7*(mJerkOld[4]+mJerk[4])));
		mPos[5] = mPosOld[5] + dtby2*((mVelOld[5]+mVel[5]) + dt7by30*((mAccOld[5]- mAcc[5]) + dtby7*(mJerkOld[5]+mJerk[5])));
		mPos[6] = mPosOld[6] + dtby2*((mVelOld[6]+mVel[6]) + dt7by30*((mAccOld[6]- mAcc[6]) + dtby7*(mJerkOld[6]+mJerk[6])));
		mPos[7] = mPosOld[7] + dtby2*((mVelOld[7]+mVel[7]) + dt7by30*((mAccOld[7]- mAcc[7]) + dtby7*(mJerkOld[7]+mJerk[7])));
		mPos[8] = mPosOld[8] + dtby2*((mVelOld[8]+mVel[8]) + dt7by30*((mAccOld[8]- mAcc[8]) + dtby7*(mJerkOld[8]+mJerk[8])));
	}
	else {
		for(int i=0; i<numArray; i++) {
			mVel[i] = mVelOld[i] + dtby2*((mAccOld[i]+mAcc[i]) + dtby6*  (mJerkOld[i]-mJerk[i]));
			mPos[i] = mPosOld[i] + dtby2*((mVelOld[i]+mVel[i]) + dt7by30*((mAccOld[i]- mAcc[i]) + dtby7*(mJerkOld[i]+mJerk[i])));
		}
	}
}
//******************************************************************
// * UpdateAccJerk function for 2 or 3 Planets 
// *(real = float for single and mixed)
// *(real = double for double)
//******************************************************************
template<class acc_real>
__device__  void UpdateAccJerk(acc_real * mPos, acc_real * mVel, acc_real* mAcc, acc_real* mJerk, int nBodies,const float * d_mass) 
{
	acc_real dx[]={0,0,0}; 
	acc_real dv[]={0,0,0}; 
	acc_real dx_back[]={0,0,0}; 
	acc_real dv_back[]={0,0,0}; 

	acc_real r2=0;
	acc_real rv=0;
	acc_real rinv=0;
	acc_real rinv3=0;

	acc_real ai0[]={0,0,0};
	acc_real ai1[]={0,0,0};
	acc_real ai2[]={0,0,0};
	acc_real ji0[]={0,0,0};
	acc_real ji1[]={0,0,0};
	acc_real ji2[]={0,0,0};

	//! planet1 and planet2
	dx[0] = mPos[6] - mPos[3]; dx_back[0] = -dx[0];
	dx[1] = mPos[7] - mPos[4]; dx_back[1] = -dx[1];
	dx[2] = mPos[8] - mPos[5]; dx_back[2] = -dx[2];
	dv[0] = mVel[6] - mVel[3]; dv_back[0] = -dv[0];
	dv[1] = mVel[7] - mVel[4]; dv_back[1] = -dv[1];
	dv[2] = mVel[8] - mVel[5]; dv_back[2] = -dv[2];

	r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
	rv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
	rinv = RSQRT(r2);
	rv *= 3.0f/r2;
	rinv *= d_mass[2];
	rinv3 =  rinv/r2;

	dx[0] *= rinv3; dx[1] *= rinv3; dx[2] *= rinv3;
	ai1[0] += dx[0]; ai1[1] += dx[1]; ai1[2] += dx[2];
	dv[0] *= rinv3; dv[1] *= rinv3; dv[2] *= rinv3;
	ji1[0] += dv[0]; ji1[1] += dv[1]; ji1[2] += dv[2];
	dx[0] *= rv; dx[1] *= rv; dx[2] *= rv;
	ji1[0] -= dx[0]; ji1[1] -= dx[1]; ji1[2] -= dx[2];


	rinv3 = rinv3/d_mass[2] * d_mass[1];

	dx_back[0] *= rinv3; dx_back[1] *= rinv3; dx_back[2] *= rinv3;
	ai2[0] += dx_back[0]; ai2[1] += dx_back[1]; ai2[2] += dx_back[2];
	dv_back[0] *= rinv3; dv_back[1] *= rinv3; dv_back[2] *= rinv3;
	ji2[0] += dv_back[0]; ji2[1] += dv_back[1]; ji2[2] += dv_back[2];
	ji2[0] -= dx_back[0]*rv; ji2[1] -= dx_back[1]*rv; ji2[2] -= dx_back[2]*rv;


	//! Star and planet 1
	dx[0] = mPos[0] - mPos[3]; dx_back[0] = -dx[0];
	dx[1] = mPos[1] - mPos[4]; dx_back[1] = -dx[1];
	dx[2] = mPos[2] - mPos[5]; dx_back[2] = -dx[2];
	dv[0] = mVel[0] - mVel[3]; dv_back[0] = -dv[0];
	dv[1] = mVel[1] - mVel[4]; dv_back[1] = -dv[1];
	dv[2] = mVel[2] - mVel[5]; dv_back[2] = -dv[2];


	r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
	rv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
	rinv = RSQRT(r2);
	rv *= 3.0f/r2;
	rinv *= d_mass[0];
	rinv3 = rinv/r2;

	dx[0] *= rinv3; dx[1] *= rinv3; dx[2] *= rinv3;
	ai1[0] += dx[0]; ai1[1] += dx[1]; ai1[2] += dx[2];
	dv[0] *= rinv3; dv[1] *= rinv3; dv[2] *= rinv3;
	ji1[0] += dv[0]; ji1[1] += dv[1]; ji1[2] += dv[2];
	dx[0] *= rv; dx[1] *= rv; dx[2] *= rv;
	ji1[0] -= dx[0]; ji1[1] -= dx[1]; ji1[2] -= dx[2];


	rinv3=rinv3/d_mass[0]*d_mass[1];

	dx_back[0] *= rinv3; dx_back[1] *= rinv3; dx_back[2] *= rinv3;
	ai0[0] += dx_back[0]; ai0[1] += dx_back[1]; ai0[2] += dx_back[2];
	dv_back[0] *= rinv3; dv_back[1] *= rinv3; dv_back[2] *= rinv3;
	ji0[0] += dv_back[0]; ji0[1] += dv_back[1]; ji0[2] += dv_back[2];
	dx_back[0] *= rv; dx_back[1] *= rv; dx_back[2] *= rv;
	ji0[0] -= dx_back[0]; ji0[1] -= dx_back[1]; ji0[2] -= dx_back[2];


	//! Star and planet 2
	dx[0] = mPos[6] - mPos[0]; dx_back[0] = -dx[0];
	dx[1] = mPos[7] - mPos[1]; dx_back[1] = -dx[1];
	dx[2] = mPos[8] - mPos[2]; dx_back[2] = -dx[2];
	dv[0] = mVel[6] - mVel[0]; dv_back[0] = -dv[0];
	dv[1] = mVel[7] - mVel[1]; dv_back[1] = -dv[1];
	dv[2] = mVel[8] - mVel[2]; dv_back[2] = -dv[2];

	r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
	rv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
	rinv = RSQRT(r2);
	rv *= 3.0f/r2;
	rinv *= d_mass[2];
	rinv3 = rinv/r2;

	dx[0] *= rinv3; dx[1] *= rinv3; dx[2] *= rinv3;
	ai0[0] += dx[0]; ai0[1] += dx[1]; ai0[2] += dx[2];
	dv[0] *= rinv3; dv[1] *= rinv3; dv[2] *= rinv3;
	ji0[0] += dv[0]; ji0[1] += dv[1]; ji0[2] += dv[2];
	dx[0] *= rv; dx[1] *= rv; dx[2] *= rv;
	ji0[0] -= dx[0]; ji0[1] -= dx[1]; ji0[2] -= dx[2];

	rinv3 = rinv3/d_mass[2] * d_mass[0];

	dx_back[0] *= rinv3; dx_back[1] *= rinv3; dx_back[2] *= rinv3;
	ai2[0] += dx_back[0]; ai2[1] += dx_back[1]; ai2[2] += dx_back[2];
	dv_back[0] *= rinv3; dv_back[1] *= rinv3; dv_back[2] *= rinv3;
	ji2[0] += dv_back[0]; ji2[1] += dv_back[1]; ji2[2] += dv_back[2];
	dx_back[0] *= rv; dx_back[1] *= rv; dx_back[2] *= rv;
	ji2[0] -= dx_back[0]; ji2[1] -= dx_back[1]; ji2[2] -= dx_back[2];

	mAcc[0] = ai0[0]; mAcc[1] = ai0[1]; mAcc[2] = ai0[2]; 
	mJerk[0] = ji0[0]; mJerk[1] = ji0[1]; mJerk[2] = ji0[2];
	mAcc[3] = ai1[0]; mAcc[4] = ai1[1]; mAcc[5] = ai1[2]; 
	mJerk[3] = ji1[0]; mJerk[4] = ji1[1]; mJerk[5] = ji1[2];
	mAcc[6] = ai2[0]; mAcc[7] = ai2[1]; mAcc[8] = ai2[2]; 
	mJerk[6] = ji2[0]; mJerk[7] = ji2[1]; mJerk[8] = ji2[2];
}

//******************************************************************
// * UpdateAccJerk function for more than 3 Planets 
// *(real = float for single and mixed)
// *(real = double for double)
// ******************************************************************/
__device__  void UpdateAccJerk_General(real * mPos, real * mVel, real* mAcc, real* mJerk, int nBodies,const float * d_mass) 
{

	real dx[]={0,0,0}; 
	real dv[]={0,0,0}; 

	{ // First calculate acceleration and jerk for the Sun
		unsigned int i = 0;
		real xi[]={mPos[i*3], mPos[i*3+1], mPos[i*3+2]};
		real vi[]={mVel[i*3], mVel[i*3+1], mVel[i*3+2]};
		real ai[]={0,0,0};
		real ji[]={0,0,0};

#pragma unroll
		for(unsigned int j=1;j<nBodies;++j)
		{
			unsigned int jj= 3*j;
			dx[0] = mPos[jj] - xi[0];
			dv[0] = mVel[jj] - vi[0];
			++jj;
			dx[1] = mPos[jj] - xi[1];
			dv[1] = mVel[jj] - vi[1];
			++jj;
			dx[2] = mPos[jj] - xi[2];
			dv[2] = mVel[jj] - vi[2];
			real r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
			real rv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
			real rinv = RSQRT(r2);
			rv *= 3./r2;
			rinv *= d_mass[j];
			real rinv3 = rinv/r2;

			//dx *= rinv3;
			//dx = rinv3*dx;
			dx[0] *= rinv3;
			dx[1] *= rinv3;
			dx[2] *= rinv3;
			//ai += dx;
			ai[0] += dx[0];
			ai[1] += dx[1];
			ai[2] += dx[2];
			//dv *= rinv3;
			//dv = rinv3*dv;
			dv[0] *= rinv3;
			dv[1] *= rinv3;
			dv[2] *= rinv3;
			//ji += dv;
			ji[0] += dv[0];
			ji[1] += dv[1];
			ji[2] += dv[2];
			//dx *= rv;
			//dx = rv*dx;
			dx[0] *= rv;
			dx[1] *= rv;
			dx[2] *= rv;
			//ji -= dx;
			ji[0] -= dx[0];
			ji[1] -= dx[1];
			ji[2] -= dx[2];
		}
		//mAcc[i] = ai;
		mAcc[i*3  ] = ai[0];
		mAcc[i*3+1] = ai[1];
		mAcc[i*3+2] = ai[2];
		//mJerk[i] = ji;
		mJerk[i*3  ] = ji[0];
		mJerk[i*3+1] = ji[1];
		mJerk[i*3+2] = ji[2];
		unsigned int ii = i*3;
		mAcc[ii  ] = ai[0];
		mJerk[ii ] = ji[0];
		++ii;
		mAcc[ii ] = ai[1];
		mJerk[ii] = ji[1];
		++ii;
		mAcc[ii ] = ai[2];
		mJerk[ii] = ji[2];
	}

#pragma unroll 
	for(unsigned int i=1;i<nBodies;++i)
	{
		//float3 xi=mPos[i];
		real xi[]={mPos[i*3], mPos[i*3+1], mPos[i*3+2]};
		real vi[]={mVel[i*3], mVel[i*3+1], mVel[i*3+2]};
		real ai[]={0,0,0};
		real ji[]={0,0,0};

#pragma unroll
		for(unsigned int j=1;j<nBodies;++j)
		{
			if(j==i) continue; // Ignore body interacting with itself
			unsigned int jj= 3*j;
			dx[0] = mPos[jj] - xi[0]; dv[0] = mVel[jj] - vi[0]; ++jj;
			dx[1] = mPos[jj] - xi[1]; dv[1] = mVel[jj] - vi[1]; ++jj;
			dx[2] = mPos[jj] - xi[2]; dv[2] = mVel[jj] - vi[2];
			//	    dx = mPos[j] - mPos[i];
			//dx = mPos[j] - xi;
			//	    dv = mVel[j] - mVel[i];
			//dv = mVel[j] - vi;
			//float r2 = dx.MagnitudeSquared();
			real r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
			//float r2 = dot(dx,dx);
			real rv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
			//float rv = dot(dx,dv);
			real rinv = RSQRT(r2);
			rv *= 3./r2;
			rinv *= d_mass[j];
			real rinv3 = rinv/r2;

			//dx *= rinv3;
			//dx = rinv3*dx;
			dx[0] = rinv3*dx[0]; dx[1] = rinv3*dx[1]; dx[2] = rinv3*dx[2];
			//ai += dx;
			ai[0] = ai[0] +dx[0]; ai[1] = ai[1] +dx[1]; ai[2] = ai[2] +dx[2];
			//dv *= rinv3;
			//dv = rinv3*dv;
			dv[0] = rinv3*dv[0]; dv[1] = rinv3*dv[1]; dv[2] = rinv3*dv[2];
			//ji += dv;
			ji[0] =ji[0] + dv[0]; ji[1] =ji[1] + dv[1]; ji[2] =ji[2] + dv[2];
			//dx *= rv;
			//dx = rv*dx;
			dx[0] = rv*dx[0]; dx[1] = rv*dx[1]; dx[2] = rv*dx[2];
			//ji -= dx;
			ji[0] = ji[0] - dx[0]; ji[1] = ji[1] - dx[1]; ji[2] = ji[2] - dx[2];
		}
		{  // But add sun's contribution last to minimize round-off error
			//	    dx = mPos[j] - mPos[i];
			dx[0] = mPos[0] - xi[0]; dx[1] = mPos[1] - xi[1]; dx[2] = mPos[2] - xi[2];
			//	    dv = mVel[j] - mVel[i];
			//dv = mVel[j] - vi;
			dv[0] = mVel[0] - vi[0]; dv[1] = mVel[1] - vi[1]; dv[2] = mVel[2] - vi[2];
			//float r2 = dx.MagnitudeSquared();
			real r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
			//float r2 = dot(dx,dx);
			real rv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
			//float rv = dot(dx,dv);
			real rinv = RSQRT(r2);
			rv *= 3./r2;
			rinv *= d_mass[0];
			real rinv3 = rinv/r2;

			//dx *= rinv3;
			//dx = rinv3*dx;
			dx[0] = rinv3*dx[0]; dx[1] = rinv3*dx[1]; dx[2] = rinv3*dx[2];
			//ai += dx;
			ai[0] = ai[0] +dx[0]; ai[1] = ai[1] +dx[1]; ai[2] = ai[2] +dx[2];
			//dv *= rinv3;
			//dv = rinv3*dv;
			dv[0] = rinv3*dv[0]; dv[1] = rinv3*dv[1]; dv[2] = rinv3*dv[2];
			//ji += dv;
			ji[0] =ji[0] + dv[0]; ji[1] =ji[1] + dv[1]; ji[2] =ji[2] + dv[2];
			//dx *= rv;
			//dx = rv*dx;
			dx[0] = rv*dx[0]; dx[1] = rv*dx[1]; dx[2] = rv*dx[2];
			//ji -= dx;
			ji[0] = ji[0] - dx[0]; ji[1] = ji[1] - dx[1]; ji[2] = ji[2] - dx[2];
		}
		unsigned int ii = i*3;
		mAcc[ii] = ai[0]; mJerk[ii] = ji[0]; ++ii;
		mAcc[ii] = ai[1]; mJerk[ii] = ji[1]; ++ii;
		mAcc[ii] = ai[2]; mJerk[ii] = ji[2];
	}
}

//__constant__ gpu_hermite_integrator_data pars;
__constant__ ensemble gpu_hermite_ens;

template<int pre>
__global__ void gpu_hermite_integrator_kernel(double dT, double h)
{
	using namespace gpu_hermite_aux;

	ensemble &ens = gpu_hermite_ens;
	int sys = threadId();
	if(sys >= ens.nsys()) { return; }

	double    T = ens.time(sys);
	double Tend = T + dT;

	//const unsigned int current_id = sys; 

	real mPos       [9];
	real mVel       [9];
	real mAcc       [9];
	real mJerk      [9];
	real mPosOld    [9];
	real mVelOld    [9];
	real mAccOld    [9];
	real mJerkOld   [9];

	float sPos       [9];
	float sVel       [9];
	float sAcc       [9];
	float sJerk      [9];
	float sPosOld    [9];
	float sVelOld    [9];
	float sAccOld    [9];
	float sJerkOld   [9];

	//const float s_mass[]={d_mass[t_start], d_mass[t_start+1],d_mass[t_start+2]};
	const float s_mass[]={ens.mass(sys, 0), ens.mass(sys,1), ens.mass(sys,2)};

	const real dtby2=h/2.;
	const real dtby3=h/3.;
	const real dtby6=h/6.;
	const real dt7by30=h*7./30.;
	const real dtby7=h*7.;
	const unsigned int nData=9;

	//load data from global memory
	mPos[0]=ens.x(sys,0);
	mPos[1]=ens.y(sys,0);
	mPos[2]=ens.z(sys,0);
	mPos[3]=ens.x(sys,1);
	mPos[4]=ens.y(sys,1);
	mPos[5]=ens.z(sys,1);
	mPos[6]=ens.x(sys,2);
	mPos[7]=ens.y(sys,2);
	mPos[8]=ens.z(sys,2);

	mVel[0]=ens.vx(sys,0);
	mVel[1]=ens.vy(sys,0);
	mVel[2]=ens.vz(sys,0);
	mVel[3]=ens.vx(sys,1);
	mVel[4]=ens.vy(sys,1);
	mVel[5]=ens.vz(sys,1);
	mVel[6]=ens.vx(sys,2);
	mVel[7]=ens.vy(sys,2);
	mVel[8]=ens.vz(sys,2);

	//UpdateAccJerk_General(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], 3, &s_mass[0]);
	if(pre==1)
		UpdateAccJerk<double>(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], 3, &s_mass[0]);
	else
	{
		doubleTofloat(sPos, mPos,9);
		doubleTofloat(sVel, mVel,9);
		UpdateAccJerk<float>(&sPos[0], &sVel[0], &sAcc[0], &sJerk[0], 3, &s_mass[0]);
		floatTodouble(mAcc,sAcc,9);
		floatTodouble(mJerk,sJerk,9);
	}

	while(T<Tend)
	{
		////Evolve(DeltaT);
		//CopyToOld();

		copyArray(mPosOld,mPos,9);
		copyArray(mVelOld,mVel,9);
		copyArray(mAccOld,mAcc,9);
		copyArray(mJerkOld,mJerk,9);
		//for(unsigned int i=0; i<nData; ++i) {
		//	mPosOld[i]=mPos[i];
		//	mVelOld[i]=mVel[i];
		//	mAccOld[i]=mAcc[i];
		//	mJerkOld[i]=mJerk[i];
		//}

		predict(mPos,mVel,mAcc,mJerk, dtby2, dtby3, h, 9);
		//for(unsigned int i=0; i<nData; ++i) {
		//	mPos[i] += h* (mVel[i]+ dtby2*(mAcc[i]+dtby3*mJerk[i]));
		//	mVel[i] += h* ( mAcc[i]+ dtby2*mJerk[i]);
		//}

		//UpdateAccJerk_General(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], 3, &s_mass[0]);
		if(pre==1)
			UpdateAccJerk<double>(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], 3, &s_mass[0]);
		else
		{
			doubleTofloat(sPos, mPos,9);
			doubleTofloat(sVel, mVel,9);
			UpdateAccJerk<float>(&sPos[0], &sVel[0], &sAcc[0], &sJerk[0], 3, &s_mass[0]);
			floatTodouble(mAcc,sAcc,9);
			floatTodouble(mJerk,sJerk,9);
		}

		//Correct(dt);
		correct(mPos,mVel,mAcc,mJerk, mPosOld,mVelOld,mAccOld,mJerkOld, dtby2, dtby6, dtby7, dt7by30, 9);
		
		//for(unsigned int i=0; i<nData; ++i) {
		//	mVel[i] = mVelOld[i] + dtby2*((mAccOld[i]+ mAcc[i]) + dtby6*  (mJerkOld[i]-mJerk[i]));
		//	mPos[i] = mPosOld[i] + dtby2*((mVelOld[i]+mVel[i]) + dt7by30*((mAccOld[i]- mAcc[i]) + dtby7*(mJerkOld[i]+mJerk[i])));
		//}

		//UpdateAccJerk_General(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], 3, &s_mass[0]);
		if(pre==1)
			UpdateAccJerk<double>(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], 3, &s_mass[0]);
		else
		{
			doubleTofloat(sPos, mPos,9);
			doubleTofloat(sVel, mVel,9);
			UpdateAccJerk<float>(&sPos[0], &sVel[0], &sAcc[0], &sJerk[0], 3, &s_mass[0]);
			floatTodouble(mAcc,sAcc,9);
			floatTodouble(mJerk,sJerk,9);
		}

		//Correct(dt);
		correct(mPos,mVel,mAcc,mJerk, mPosOld,mVelOld,mAccOld,mJerkOld, dtby2, dtby6, dtby7, dt7by30, 9);
		//for(unsigned int i=0; i<nData; ++i) {
		//	mVel[i] = mVelOld[i] + dtby2*((mAccOld[i]+ mAcc[i]) + dtby6*  (mJerkOld[i]-mJerk[i]));
		//	mPos[i] = mPosOld[i] + dtby2*((mVelOld[i]+mVel[i]) + dt7by30*((mAccOld[i]- mAcc[i]) + dtby7*(mJerkOld[i]+mJerk[i])));
		//}

		T += h;
	}
	ens.x(sys,0)=mPos[0];
	ens.y(sys,0)=mPos[1];
	ens.z(sys,0)=mPos[2];
	ens.x(sys,1)=mPos[3];
	ens.y(sys,1)=mPos[4];
	ens.z(sys,1)=mPos[5];
	ens.x(sys,2)=mPos[6];
	ens.y(sys,2)=mPos[7];
	ens.z(sys,2)=mPos[8];

	ens.vx(sys,0)=mVel[0];
	ens.vy(sys,0)=mVel[1];
	ens.vz(sys,0)=mVel[2];
	ens.vx(sys,1)=mVel[3];
	ens.vy(sys,1)=mVel[4];
	ens.vz(sys,1)=mVel[5];
	ens.vx(sys,2)=mVel[6];
	ens.vy(sys,2)=mVel[7];
	ens.vz(sys,2)=mVel[8];
}

void gpu_hermite_integrator::integrate(gpu_ensemble &ens, double dT)
{
	// Upload the kernel parameters
	if(ens.last_integrator() != this)
	{
		ens.set_last_integrator(this);
		configure_grid(gridDim, threadsPerBlock, ens.nsys());

		cudaMemcpyToSymbol(gpu_hermite_ens, &ens, sizeof(gpu_hermite_ens));
	}

	// execute the kernel
	switch(prec){
		// double precision
		case 1:
			gpu_hermite_integrator_kernel<1><<<gridDim, threadsPerBlock>>>(dT, h);
			break;
		// signle precision
		case 2:
			gpu_hermite_integrator_kernel<2><<<gridDim, threadsPerBlock>>>(dT, h);
			break;
		// mixed precision
		case 3:
			gpu_hermite_integrator_kernel<3><<<gridDim, threadsPerBlock>>>(dT, h);
			break;
	}

}

