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

#define XYZ 0
//#define MAX_STEPS 10000
//#if PRECISION == 3
// Double precision
#define RSQRT(x) rsqrt(x)
#define SQRT(x)   sqrt(x)
//typedef double real;
//#else
//// Single or Mixed precision
//#define RSQRT(x) rsqrtf(x)
//#define SQRT(x)   sqrtf(x)
//typedef float real;
//#endif


//******************************************************************
// * UpdateAccJerk function for 2 or 3 Planets 
// *(real = float for single and mixed)
// *(real = double for double)
//******************************************************************
__device__  void UpdateAccJerk(real * mPos, real * mVel, real* mAcc, real* mJerk, int nBodies,const float * d_mass) 
{
	real dx[]={0,0,0}; 
	real dv[]={0,0,0}; 
	real dx_back[]={0,0,0}; 
	real dv_back[]={0,0,0}; 

	real r2=0;
	real rv=0;
	real rinv=0;
	real rinv3=0;

	real ai0[]={0,0,0};
	real ai1[]={0,0,0};
	real ai2[]={0,0,0};
	real ji0[]={0,0,0};
	real ji1[]={0,0,0};
	real ji2[]={0,0,0};

	//! planet1 and planet2
#if XYZ
	dx[0] = mPos[6] - mPos[3]; dx_back[0] = -dx[0];
	dx[1] = mPos[7] - mPos[4]; dx_back[1] = -dx[1];
	dx[2] = mPos[8] - mPos[5]; dx_back[2] = -dx[2];
	dv[0] = mVel[6] - mVel[3]; dv_back[0] = -dv[0];
	dv[1] = mVel[7] - mVel[4]; dv_back[1] = -dv[1];
	dv[2] = mVel[8] - mVel[5]; dv_back[2] = -dv[2];
#else
	dx[0] = mPos[2] - mPos[1]; dx_back[0] = -dx[0];
	dx[1] = mPos[5] - mPos[4]; dx_back[1] = -dx[1];
	dx[2] = mPos[8] - mPos[7]; dx_back[2] = -dx[2];
	dv[0] = mVel[2] - mVel[1]; dv_back[0] = -dv[0];
	dv[1] = mVel[5] - mVel[4]; dv_back[1] = -dv[1];
	dv[2] = mVel[8] - mVel[7]; dv_back[2] = -dv[2];
#endif

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
#if XYZ
	dx[0] = mPos[0] - mPos[3]; dx_back[0] = -dx[0];
	dx[1] = mPos[1] - mPos[4]; dx_back[1] = -dx[1];
	dx[2] = mPos[2] - mPos[5]; dx_back[2] = -dx[2];
	dv[0] = mVel[0] - mVel[3]; dv_back[0] = -dv[0];
	dv[1] = mVel[1] - mVel[4]; dv_back[1] = -dv[1];
	dv[2] = mVel[2] - mVel[5]; dv_back[2] = -dv[2];
#else
	dx[0] = mPos[0] - mPos[1]; dx_back[0] = -dx[0];
	dx[1] = mPos[3] - mPos[4]; dx_back[1] = -dx[1];
	dx[2] = mPos[6] - mPos[7]; dx_back[2] = -dx[2];
	dv[0] = mVel[0] - mVel[1]; dv_back[0] = -dv[0];
	dv[1] = mVel[3] - mVel[4]; dv_back[1] = -dv[1];
	dv[2] = mVel[6] - mVel[7]; dv_back[2] = -dv[2];
#endif


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
#if XYZ
	dx[0] = mPos[6] - mPos[0]; dx_back[0] = -dx[0];
	dx[1] = mPos[7] - mPos[1]; dx_back[1] = -dx[1];
	dx[2] = mPos[8] - mPos[2]; dx_back[2] = -dx[2];
	dv[0] = mVel[6] - mVel[0]; dv_back[0] = -dv[0];
	dv[1] = mVel[7] - mVel[1]; dv_back[1] = -dv[1];
	dv[2] = mVel[8] - mVel[2]; dv_back[2] = -dv[2];
#else
	dx[0] = mPos[2] - mPos[0]; dx_back[0] = -dx[0];
	dx[1] = mPos[5] - mPos[3]; dx_back[1] = -dx[1];
	dx[2] = mPos[8] - mPos[6]; dx_back[2] = -dx[2];
	dv[0] = mVel[2] - mVel[0]; dv_back[0] = -dv[0];
	dv[1] = mVel[5] - mVel[3]; dv_back[1] = -dv[1];
	dv[2] = mVel[8] - mVel[6]; dv_back[2] = -dv[2];
#endif

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

#if XYZ
	mAcc[0] = ai0[0]; mAcc[1] = ai0[1]; mAcc[2] = ai0[2]; 
	mJerk[0] = ji0[0]; mJerk[1] = ji0[1]; mJerk[2] = ji0[2];
	mAcc[3] = ai1[0]; mAcc[4] = ai1[1]; mAcc[5] = ai1[2]; 
	mJerk[3] = ji1[0]; mJerk[4] = ji1[1]; mJerk[5] = ji1[2];
	mAcc[6] = ai2[0]; mAcc[7] = ai2[1]; mAcc[8] = ai2[2]; 
	mJerk[6] = ji2[0]; mJerk[7] = ji2[1]; mJerk[8] = ji2[2];
#else
	mAcc[0] = ai0[0]; mAcc[3] = ai0[1]; mAcc[6] = ai0[2]; 
	mJerk[0] = ji0[0]; mJerk[3] = ji0[1]; mJerk[6] = ji0[2];
	mAcc[1] = ai1[0]; mAcc[4] = ai1[1]; mAcc[7] = ai1[2]; 
	mJerk[1] = ji1[0]; mJerk[4] = ji1[1]; mJerk[7] = ji1[2];
	mAcc[2] = ai2[0]; mAcc[5] = ai2[1]; mAcc[8] = ai2[2]; 
	mJerk[2] = ji2[0]; mJerk[5] = ji2[1]; mJerk[8] = ji2[2];
#endif
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

__global__ void gpu_hermite_integrator_kernel(double dT, double h)
{
	using namespace gpu_hermite_aux;

	ensemble &ens = gpu_hermite_ens;
	int sys = threadId();
	if(sys >= ens.nsys()) { return; }

	double    T = ens.time(sys);
	double Tend = T + dT;

	//const unsigned int current_id = sys; 

	real mPos       [(2+1)*3];
	real mVel       [(2+1)*3];
	real mAcc       [(2+1)*3];
	real mJerk      [(2+1)*3];
	real mPosOld    [(2+1)*3];
	real mVelOld    [(2+1)*3];
	real mAccOld    [(2+1)*3];
	real mJerkOld   [(2+1)*3];

	//const float s_mass[]={d_mass[t_start], d_mass[t_start+1],d_mass[t_start+2]};
	const float s_mass[]={ens.mass(sys, 0), ens.mass(sys,1), ens.mass(sys,2)};

	const real dtby2=h/2.;
	const real dtby3=h/3.;
	const real dtby6=h/6.;
	const real dt7by30=h*7./30.;
	const real dtby7=h*7.;

	//load data from global memory

	//const unsigned int nData=nBodies*3;
	const unsigned int nData=3*3;
#if XYZ
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
#else
	mPos[0]=ens.x(sys,0);
	mPos[1]=ens.x(sys,1);
	mPos[2]=ens.x(sys,2);
	mPos[3]=ens.y(sys,0);
	mPos[4]=ens.y(sys,1);
	mPos[5]=ens.y(sys,2);
	mPos[6]=ens.z(sys,0);
	mPos[7]=ens.z(sys,1);
	mPos[8]=ens.z(sys,2);

	mVel[0]=ens.vx(sys,0);
	mVel[1]=ens.vx(sys,1);
	mVel[2]=ens.vx(sys,2);
	mVel[3]=ens.vy(sys,0);
	mVel[4]=ens.vy(sys,1);
	mVel[5]=ens.vy(sys,2);
	mVel[6]=ens.vz(sys,0);
	mVel[7]=ens.vz(sys,1);
	mVel[8]=ens.vz(sys,2);
#endif

	//UpdateAccJerk_General(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], 3, &s_mass[0]);
	UpdateAccJerk(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], 3, &s_mass[0]);

	while(T<Tend)
	{
		////Evolve(DeltaT);
		//CopyToOld();

		for(unsigned int i=0; i<nData; ++i) {
			mPosOld[i]=mPos[i];
			mVelOld[i]=mVel[i];
			mAccOld[i]=mAcc[i];
			mJerkOld[i]=mJerk[i];
		}

		for(unsigned int i=0; i<nData; ++i) {
			mPos[i] += h* (mVel[i]+ dtby2*(mAcc[i]+dtby3*mJerk[i]));
			mVel[i] += h* ( mAcc[i]+ dtby2*mJerk[i]);
		}

		//UpdateAccJerk_General(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], 3, &s_mass[0]);
		UpdateAccJerk(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], 3, &s_mass[0]);
                
		//Correct(dt);
		for(unsigned int i=0; i<nData; ++i) {
			mVel[i] = mVelOld[i] + dtby2*((mAccOld[i]+ mAcc[i]) + dtby6*  (mJerkOld[i]-mJerk[i]));
			mPos[i] = mPosOld[i] + dtby2*((mVelOld[i]+mVel[i]) + dt7by30*((mAccOld[i]- mAcc[i]) + dtby7*(mJerkOld[i]+mJerk[i])));
		}

		//UpdateAccJerk_General(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], 3, &s_mass[0]);
		UpdateAccJerk(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], 3, &s_mass[0]);

		//Correct(dt);
		for(unsigned int i=0; i<nData; ++i) {
			mVel[i] = mVelOld[i] + dtby2*((mAccOld[i]+ mAcc[i]) + dtby6*  (mJerkOld[i]-mJerk[i]));
			mPos[i] = mPosOld[i] + dtby2*((mVelOld[i]+mVel[i]) + dt7by30*((mAccOld[i]- mAcc[i]) + dtby7*(mJerkOld[i]+mJerk[i])));
		}

		T += h;
	}
#if XYZ

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
#else

	ens.x(sys,0)=mPos[0];
	ens.x(sys,1)=mPos[1];
	ens.x(sys,2)=mPos[2];
	ens.y(sys,0)=mPos[3];
	ens.y(sys,1)=mPos[4];
	ens.y(sys,2)=mPos[5];
	ens.z(sys,0)=mPos[6];
	ens.z(sys,1)=mPos[7];
	ens.z(sys,2)=mPos[8];

	ens.vx(sys,0)=mVel[0];
	ens.vx(sys,1)=mVel[1];
	ens.vx(sys,2)=mVel[2];
	ens.vy(sys,0)=mVel[3];
	ens.vy(sys,1)=mVel[4];
	ens.vy(sys,2)=mVel[5];
	ens.vz(sys,0)=mVel[6];
	ens.vz(sys,1)=mVel[7];
	ens.vz(sys,2)=mVel[8];
        ens.time(sys)=T;

#endif
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
	gpu_hermite_integrator_kernel<<<gridDim, threadsPerBlock>>>(dT, h);
}

