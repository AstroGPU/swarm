#include "swarm.h"
#include "hermite_gpu.h"

namespace swarm {
namespace hermite_gpu {

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

	// Specify position/velocity type to be used for various modes
	template<int pre>
	struct pos_type {}; 
	template<>
	struct pos_type<1> { typedef double type; };
	template<>
	struct pos_type<2> { typedef float type; };
	template<>
	struct pos_type<3> { typedef double type; };

	// Specify acceleration/jerk type to be used for various modes
	template<int pre>
	struct acc_type {}; 
	template<>
	struct acc_type<1> { typedef double type; };
	template<>
	struct acc_type<2> { typedef float type; };
	template<>
	struct acc_type<3> { typedef float type; };

}

// Be careful, do we want the lower precission sqrt?
#define RSQRT(x) rsqrtf(x)
#define SQRT(x)   sqrtf(x)

template<unsigned int N, typename destT, typename srcT>
inline __device__ void copyArray(destT *target, srcT *source)
{
	if(N>0)  target[0]=source[0];
	if(N>1)  target[1]=source[1];
	if(N>2)  target[2]=source[2];
	if(N>3)  target[3]=source[3];
	if(N>4)  target[4]=source[4];
	if(N>5)  target[5]=source[5];
	if(N>6)  target[6]=source[6];
	if(N>7)  target[7]=source[7];
	if(N>8)  target[8]=source[8];
	if(N>9)  target[9]=source[9];
	if(N>10) target[10]=source[10];
	if(N>11) target[11]=source[11];
	if(N>12) target[12]=source[12];
	if(N>13) target[13]=source[13];
	if(N>14) target[14]=source[14];
	if(N>15) target[15]=source[15];
	if(N>16) target[16]=source[16];
	if(N>17) target[17]=source[17];
	if(N>18) {
		 for(int i=18; i<N; i++)
		    target[i]=source[i];
		 }
}


template<unsigned int N, typename real_hi, typename real_lo>
inline __device__ void predict(real_hi *mPos, real_hi *mVel, real_lo *mAcc, real_lo *mJerk, const real_lo dtby2, const real_lo dtby3, real_hi h)
{
	if(N>0)	mPos[0] += h* (mVel[0]+ dtby2*(mAcc[0]+dtby3*mJerk[0]));
	if(N>1)	mPos[1] += h* (mVel[1]+ dtby2*(mAcc[1]+dtby3*mJerk[1]));
	if(N>2)	mPos[2] += h* (mVel[2]+ dtby2*(mAcc[2]+dtby3*mJerk[2]));
	if(N>3)	mPos[3] += h* (mVel[3]+ dtby2*(mAcc[3]+dtby3*mJerk[3]));
	if(N>4)	mPos[4] += h* (mVel[4]+ dtby2*(mAcc[4]+dtby3*mJerk[4]));
	if(N>5)	mPos[5] += h* (mVel[5]+ dtby2*(mAcc[5]+dtby3*mJerk[5]));
	if(N>6)	mPos[6] += h* (mVel[6]+ dtby2*(mAcc[6]+dtby3*mJerk[6]));
	if(N>7)	mPos[7] += h* (mVel[7]+ dtby2*(mAcc[7]+dtby3*mJerk[7]));
	if(N>8)	mPos[8] += h* (mVel[8]+ dtby2*(mAcc[8]+dtby3*mJerk[8]));
	if(N>9)	mPos[9] += h* (mVel[9]+ dtby2*(mAcc[9]+dtby3*mJerk[9]));
	if(N>10) mPos[10] += h* (mVel[10]+ dtby2*(mAcc[10]+dtby3*mJerk[10]));
	if(N>11) mPos[11] += h* (mVel[11]+ dtby2*(mAcc[11]+dtby3*mJerk[11]));
	if(N>12) mPos[12] += h* (mVel[12]+ dtby2*(mAcc[12]+dtby3*mJerk[12]));
	if(N>13) mPos[13] += h* (mVel[13]+ dtby2*(mAcc[13]+dtby3*mJerk[13]));
	if(N>14) mPos[14] += h* (mVel[14]+ dtby2*(mAcc[14]+dtby3*mJerk[14]));
	if(N>15) mPos[15] += h* (mVel[15]+ dtby2*(mAcc[15]+dtby3*mJerk[15]));
	if(N>16) mPos[16] += h* (mVel[16]+ dtby2*(mAcc[16]+dtby3*mJerk[16]));
	if(N>17) mPos[17] += h* (mVel[17]+ dtby2*(mAcc[17]+dtby3*mJerk[17]));
	if(N>0)	mVel[0] += h* (mAcc[0]+ dtby2*mJerk[0]);
	if(N>1)	mVel[1] += h* (mAcc[1]+ dtby2*mJerk[1]);
	if(N>2)	mVel[2] += h* (mAcc[2]+ dtby2*mJerk[2]);
	if(N>3)	mVel[3] += h* (mAcc[3]+ dtby2*mJerk[3]);
	if(N>4)	mVel[4] += h* (mAcc[4]+ dtby2*mJerk[4]);
	if(N>5)	mVel[5] += h* (mAcc[5]+ dtby2*mJerk[5]);
	if(N>6)	mVel[6] += h* (mAcc[6]+ dtby2*mJerk[6]);
	if(N>7)	mVel[7] += h* (mAcc[7]+ dtby2*mJerk[7]);
	if(N>8)	mVel[8] += h* (mAcc[8]+ dtby2*mJerk[8]);
	if(N>9)	mVel[9] += h* (mAcc[9]+ dtby2*mJerk[9]);
	if(N>10) mVel[10] += h* (mAcc[10]+ dtby2*mJerk[10]);
	if(N>11) mVel[11] += h* (mAcc[11]+ dtby2*mJerk[11]);
	if(N>12) mVel[12] += h* (mAcc[12]+ dtby2*mJerk[12]);
	if(N>13) mVel[13] += h* (mAcc[13]+ dtby2*mJerk[13]);
	if(N>14) mVel[14] += h* (mAcc[14]+ dtby2*mJerk[14]);
	if(N>15) mVel[15] += h* (mAcc[15]+ dtby2*mJerk[15]);
	if(N>16) mVel[16] += h* (mAcc[16]+ dtby2*mJerk[16]);
	if(N>17) mVel[17] += h* (mAcc[17]+ dtby2*mJerk[17]);
	if(N>18)
		{	
		for(int i=18; i<N; i++) {
			mPos[i] += h* (mVel[i]+ dtby2*(mAcc[i]+dtby3*mJerk[i]));
			mVel[i] += h* (mAcc[i]+ dtby2*mJerk[i]);
			}
		}
}

template<unsigned int N, typename real_hi, typename real_lo>
inline __device__ void correct(real_hi *mPos, real_hi *mVel, real_lo *mAcc, real_lo *mJerk, 
		real_hi *mPosOld, real_hi *mVelOld, real_lo *mAccOld, real_lo *mJerkOld, 
		const real_lo dtby2, const real_lo dtby6, const real_lo dtby7, const real_lo dt7by30)
{
	if(N>0)	mVel[0] = mVelOld[0] + dtby2*((mAccOld[0]+mAcc[0]) + dtby6*  (mJerkOld[0]-mJerk[0]));
	if(N>1)	mVel[1] = mVelOld[1] + dtby2*((mAccOld[1]+mAcc[1]) + dtby6*  (mJerkOld[1]-mJerk[1]));
	if(N>2)	mVel[2] = mVelOld[2] + dtby2*((mAccOld[2]+mAcc[2]) + dtby6*  (mJerkOld[2]-mJerk[2]));
	if(N>3)	mVel[3] = mVelOld[3] + dtby2*((mAccOld[3]+mAcc[3]) + dtby6*  (mJerkOld[3]-mJerk[3]));
	if(N>4)	mVel[4] = mVelOld[4] + dtby2*((mAccOld[4]+mAcc[4]) + dtby6*  (mJerkOld[4]-mJerk[4]));
	if(N>5)	mVel[5] = mVelOld[5] + dtby2*((mAccOld[5]+mAcc[5]) + dtby6*  (mJerkOld[5]-mJerk[5]));
	if(N>6)	mVel[6] = mVelOld[6] + dtby2*((mAccOld[6]+mAcc[6]) + dtby6*  (mJerkOld[6]-mJerk[6]));
	if(N>7)	mVel[7] = mVelOld[7] + dtby2*((mAccOld[7]+mAcc[7]) + dtby6*  (mJerkOld[7]-mJerk[7]));
	if(N>8)	mVel[8] = mVelOld[8] + dtby2*((mAccOld[8]+mAcc[8]) + dtby6*  (mJerkOld[8]-mJerk[8]));
	if(N>9)	mVel[9] = mVelOld[9] + dtby2*((mAccOld[9]+mAcc[9]) + dtby6*  (mJerkOld[9]-mJerk[9]));
	if(N>10)	mVel[10] = mVelOld[10] + dtby2*((mAccOld[10]+mAcc[10]) + dtby6*  (mJerkOld[10]-mJerk[10]));
	if(N>11)	mVel[11] = mVelOld[11] + dtby2*((mAccOld[11]+mAcc[11]) + dtby6*  (mJerkOld[11]-mJerk[11]));
	if(N>12)	mVel[12] = mVelOld[12] + dtby2*((mAccOld[12]+mAcc[12]) + dtby6*  (mJerkOld[12]-mJerk[12]));
	if(N>13)	mVel[13] = mVelOld[13] + dtby2*((mAccOld[13]+mAcc[13]) + dtby6*  (mJerkOld[13]-mJerk[13]));
	if(N>14)	mVel[14] = mVelOld[14] + dtby2*((mAccOld[14]+mAcc[14]) + dtby6*  (mJerkOld[14]-mJerk[14]));
	if(N>15)	mVel[15] = mVelOld[15] + dtby2*((mAccOld[15]+mAcc[15]) + dtby6*  (mJerkOld[15]-mJerk[15]));
	if(N>16)	mVel[16] = mVelOld[16] + dtby2*((mAccOld[16]+mAcc[16]) + dtby6*  (mJerkOld[16]-mJerk[16]));
	if(N>17)	mVel[17] = mVelOld[17] + dtby2*((mAccOld[17]+mAcc[17]) + dtby6*  (mJerkOld[17]-mJerk[17]));
	if(N>0)	mPos[0] = mPosOld[0] + dtby2*((mVelOld[0]+mVel[0]) + dt7by30*((mAccOld[0]- mAcc[0]) + dtby7*(mJerkOld[0]+mJerk[0])));
	if(N>1)	mPos[1] = mPosOld[1] + dtby2*((mVelOld[1]+mVel[1]) + dt7by30*((mAccOld[1]- mAcc[1]) + dtby7*(mJerkOld[1]+mJerk[1])));
	if(N>2)	mPos[2] = mPosOld[2] + dtby2*((mVelOld[2]+mVel[2]) + dt7by30*((mAccOld[2]- mAcc[2]) + dtby7*(mJerkOld[2]+mJerk[2])));
	if(N>3)	mPos[3] = mPosOld[3] + dtby2*((mVelOld[3]+mVel[3]) + dt7by30*((mAccOld[3]- mAcc[3]) + dtby7*(mJerkOld[3]+mJerk[3])));
	if(N>4)	mPos[4] = mPosOld[4] + dtby2*((mVelOld[4]+mVel[4]) + dt7by30*((mAccOld[4]- mAcc[4]) + dtby7*(mJerkOld[4]+mJerk[4])));
	if(N>5)	mPos[5] = mPosOld[5] + dtby2*((mVelOld[5]+mVel[5]) + dt7by30*((mAccOld[5]- mAcc[5]) + dtby7*(mJerkOld[5]+mJerk[5])));
	if(N>6)	mPos[6] = mPosOld[6] + dtby2*((mVelOld[6]+mVel[6]) + dt7by30*((mAccOld[6]- mAcc[6]) + dtby7*(mJerkOld[6]+mJerk[6])));
	if(N>7)	mPos[7] = mPosOld[7] + dtby2*((mVelOld[7]+mVel[7]) + dt7by30*((mAccOld[7]- mAcc[7]) + dtby7*(mJerkOld[7]+mJerk[7])));
	if(N>8)	mPos[8] = mPosOld[8] + dtby2*((mVelOld[8]+mVel[8]) + dt7by30*((mAccOld[8]- mAcc[8]) + dtby7*(mJerkOld[8]+mJerk[8])));
	if(N>9)	mPos[9] = mPosOld[9] + dtby2*((mVelOld[9]+mVel[9]) + dt7by30*((mAccOld[9]- mAcc[9]) + dtby7*(mJerkOld[9]+mJerk[9])));
	if(N>10)	mPos[10] = mPosOld[10] + dtby2*((mVelOld[10]+mVel[10]) + dt7by30*((mAccOld[10]- mAcc[10]) + dtby7*(mJerkOld[10]+mJerk[10])));
	if(N>11)	mPos[11] = mPosOld[11] + dtby2*((mVelOld[11]+mVel[11]) + dt7by30*((mAccOld[11]- mAcc[11]) + dtby7*(mJerkOld[11]+mJerk[11])));
	if(N>12)	mPos[12] = mPosOld[12] + dtby2*((mVelOld[12]+mVel[12]) + dt7by30*((mAccOld[12]- mAcc[12]) + dtby7*(mJerkOld[12]+mJerk[12])));
	if(N>13)	mPos[13] = mPosOld[13] + dtby2*((mVelOld[13]+mVel[13]) + dt7by30*((mAccOld[13]- mAcc[13]) + dtby7*(mJerkOld[13]+mJerk[13])));
	if(N>14)	mPos[14] = mPosOld[14] + dtby2*((mVelOld[14]+mVel[14]) + dt7by30*((mAccOld[14]- mAcc[14]) + dtby7*(mJerkOld[14]+mJerk[14])));
	if(N>15)	mPos[15] = mPosOld[15] + dtby2*((mVelOld[15]+mVel[15]) + dt7by30*((mAccOld[15]- mAcc[15]) + dtby7*(mJerkOld[15]+mJerk[15])));
	if(N>16)	mPos[16] = mPosOld[16] + dtby2*((mVelOld[16]+mVel[16]) + dt7by30*((mAccOld[16]- mAcc[16]) + dtby7*(mJerkOld[16]+mJerk[16])));
	if(N>17)	mPos[17] = mPosOld[17] + dtby2*((mVelOld[17]+mVel[17]) + dt7by30*((mAccOld[17]- mAcc[17]) + dtby7*(mJerkOld[17]+mJerk[17])));
	if(N>18)
	    {
		for(int i=18; i<N; i++) {
			mVel[i] = mVelOld[i] + dtby2*((mAccOld[i]+mAcc[i]) + dtby6*  (mJerkOld[i]-mJerk[i]));
			mPos[i] = mPosOld[i] + dtby2*((mVelOld[i]+mVel[i]) + dt7by30*((mAccOld[i]- mAcc[i]) + dtby7*(mJerkOld[i]+mJerk[i])));
		}	
	    }


}

//******************************************************************
// * UpdateAccJerk function for 2 or 3 Planets 
// *(real_hi = double)
// *(real_lo = float for single and mixed)
// *(real_lo = double for double)
//******************************************************************
template<unsigned int nBodies, typename real_hi, typename real_lo>
__device__  void UpdateAccJerk23(real_hi * mPos, real_hi * mVel, real_lo* mAcc, real_lo* mJerk, const float * d_mass) 
{
	// Need to enforce that that nbod == 2 or 3 at compile time
	if((nBodies<3)||(nBodies>4)) 
	{ return; }

	real_hi dx[]={0,0,0}; 
	real_hi dv[]={0,0,0}; 
	real_hi dx_back[]={0,0,0}; 
	real_hi dv_back[]={0,0,0}; 

	real_hi r2=0;
	real_hi rv=0;
	real_hi rinv=0;
	real_hi rinv3=0;

	real_lo ai0[]={0,0,0};
	real_lo ai1[]={0,0,0};
	real_lo ai2[]={0,0,0};
	real_lo ji0[]={0,0,0};
	real_lo ji1[]={0,0,0};
	real_lo ji2[]={0,0,0};

	//if(nBodies ==4) {
	real_lo ai3[]={0,0,0};
	real_lo ji3[]={0,0,0};
	//}

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


	if(nBodies==4) {
		//! planet1 and planet 3
		dx[0] = mPos[9] - mPos[3]; dx_back[0] = -dx[0];
		dx[1] = mPos[10] - mPos[4]; dx_back[1] = -dx[1];
		dx[2] = mPos[11] - mPos[5]; dx_back[2] = -dx[2];
		dv[0] = mVel[9] - mVel[3]; dv_back[0] = -dv[0];
		dv[1] = mVel[10] - mVel[4]; dv_back[1] = -dv[1];
		dv[2] = mVel[11] - mVel[5]; dv_back[2] = -dv[2];

		r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
		rv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
		rinv = RSQRT(r2);
		rv *= 3.0f/r2;
		rinv *= d_mass[3];
		rinv3 = rinv/r2;

		dx[0] *= rinv3; dx[1] *= rinv3; dx[2] *= rinv3;
		ai1[0] += dx[0]; ai1[1] += dx[1]; ai1[2] += dx[2];
		dv[0] *= rinv3; dv[1] *= rinv3; dv[2] *= rinv3;
		ji1[0] += dv[0]; ji1[1] += dv[1]; ji1[2] += dv[2];
		dx[0] *= rv; dx[1] *= rv; dx[2] *= rv;
		ji1[0] -= dx[0]; ji1[1] -= dx[1]; ji1[2] -= dx[2];


		rinv3 = rinv3/d_mass[3] * d_mass[1];

		dx_back[0] *= rinv3; dx_back[1] *= rinv3; dx_back[2] *= rinv3;
		ai3[0] += dx_back[0]; ai3[1] += dx_back[1]; ai3[2] += dx_back[2];
		dv_back[0] *= rinv3; dv_back[1] *= rinv3; dv_back[2] *= rinv3;
		ji3[0] += dv_back[0]; ji3[1] += dv_back[1]; ji3[2] += dv_back[2];
		dx_back[0] *= rv; dx_back[1] *= rv; dx_back[2] *= rv;
		ji3[0] -= dx_back[0]; ji3[1] -= dx_back[1]; ji3[2] -= dx_back[2];

		//! planet2 and planet 3
		dx[0] = mPos[9] - mPos[6]; dx_back[0] = -dx[0];
		dx[1] = mPos[10] - mPos[7]; dx_back[1] = -dx[1];
		dx[2] = mPos[11] - mPos[8]; dx_back[2] = -dx[2];
		dv[0] = mVel[9] - mVel[6]; dv_back[0] = -dv[0];
		dv[1] = mVel[10] - mVel[7]; dv_back[1] = -dv[1];
		dv[2] = mVel[11] - mVel[8]; dv_back[2] = -dv[2];


		r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
		rv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
		rinv = RSQRT(r2);
		rv *= 3.0f/r2;
		rinv *= d_mass[3];
		rinv3 = rinv/r2;

		dx[0] *= rinv3; dx[1] *= rinv3; dx[2] *= rinv3;
		ai2[0] += dx[0]; ai2[1] += dx[1]; ai2[2] += dx[2];
		dv[0] *= rinv3; dv[1] *= rinv3; dv[2] *= rinv3;
		ji2[0] += dv[0]; ji2[1] += dv[1]; ji2[2] += dv[2];
		dx[0] *= rv; dx[1] *= rv; dx[2] *= rv;
		ji2[0] -= dx[0]; ji2[1] -= dx[1]; ji2[2] -= dx[2];

		rinv3 = rinv3/d_mass[3] * d_mass[2];

		dx_back[0] *= rinv3; dx_back[1] *= rinv3; dx_back[2] *= rinv3;
		ai3[0] += dx_back[0]; ai3[1] += dx_back[1]; ai3[2] += dx_back[2];
		dv_back[0] *= rinv3; dv_back[1] *= rinv3; dv_back[2] *= rinv3;
		ji3[0] += dv_back[0]; ji3[1] += dv_back[1]; ji3[2] += dv_back[2];
		dx_back[0] *= rv; dx_back[1] *= rv; dx_back[2] *= rv;
		ji3[0] -= dx_back[0]; ji3[1] -= dx_back[1]; ji3[2] -= dx_back[2];
	}


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

	
	if(nBodies==4){
		//! Star and planet 3
		dx[0] = mPos[9] - mPos[0]; dx_back[0] = -dx[0];
		dx[1] = mPos[10] - mPos[1]; dx_back[1] = -dx[1];
		dx[2] = mPos[11] - mPos[2]; dx_back[2] = -dx[2];
		dv[0] = mVel[9] - mVel[0]; dv_back[0] = -dv[0];
		dv[1] = mVel[10] - mVel[1]; dv_back[1] = -dv[1];
		dv[2] = mVel[11] - mVel[2]; dv_back[2] = -dv[2];

		r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
		rv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
		rinv = RSQRT(r2);
		rv *= 3.0f/r2;
		rinv *= d_mass[3];
		rinv3 = rinv/r2;

		dx[0] *= rinv3; dx[1] *= rinv3; dx[2] *= rinv3;
		ai0[0] += dx[0]; ai0[1] += dx[1]; ai0[2] += dx[2];
		dv[0] *= rinv3; dv[1] *= rinv3; dv[2] *= rinv3;
		ji0[0] += dv[0]; ji0[1] += dv[1]; ji0[2] += dv[2];
		dx[0] *= rv; dx[1] *= rv; dx[2] *= rv;
		ji0[0] -= dx[0]; ji0[1] -= dx[1]; ji0[2] -= dx[2];

		rinv3 = rinv3/d_mass[3] * d_mass[0];

		dx_back[0] *= rinv3; dx_back[1] *= rinv3; dx_back[2] *= rinv3;
		ai3[0] += dx_back[0]; ai3[1] += dx_back[1]; ai3[2] += dx_back[2];
		dv_back[0] *= rinv3; dv_back[1] *= rinv3; dv_back[2] *= rinv3;
		ji3[0] += dv_back[0]; ji3[1] += dv_back[1]; ji3[2] += dv_back[2];
		dx_back[0] *= rv; dx_back[1] *= rv; dx_back[2] *= rv;
		ji3[0] -= dx_back[0]; ji3[1] -= dx_back[1]; ji3[2] -= dx_back[2];

		mAcc[9] = ai3[0]; mAcc[10] = ai3[1]; mAcc[11] = ai3[2]; 
		mJerk[9] = ji3[0]; mJerk[10] = ji3[1]; mJerk[11] = ji3[2];
	}
	mAcc[0] = ai0[0]; mAcc[1] = ai0[1]; mAcc[2] = ai0[2]; 
	mJerk[0] = ji0[0]; mJerk[1] = ji0[1]; mJerk[2] = ji0[2];
	mAcc[3] = ai1[0]; mAcc[4] = ai1[1]; mAcc[5] = ai1[2]; 
	mJerk[3] = ji1[0]; mJerk[4] = ji1[1]; mJerk[5] = ji1[2];
	mAcc[6] = ai2[0]; mAcc[7] = ai2[1]; mAcc[8] = ai2[2]; 
	mJerk[6] = ji2[0]; mJerk[7] = ji2[1]; mJerk[8] = ji2[2];
}

//******************************************************************
// * UpdateAccJerk function for more than 3 Planets 
// *(real_hi = double)
// *(real_lo = float for single and mixed)
// *(real_lo = double for double)
// ******************************************************************/
template<unsigned int nBodies, typename real_hi, typename real_lo>
//template<unsigned int nBodies>
__device__  void UpdateAccJerkGeneral(real_hi * mPos, real_hi * mVel, real_lo* mAcc, real_lo* mJerk, const float * d_mass) 
{

	real_hi dx[]={0,0,0}; 
	real_hi dv[]={0,0,0}; 

	{ // First calculate acceleration and jerk for the Sun
		unsigned int i = 0;
		real_hi xi[]={mPos[i*3], mPos[i*3+1], mPos[i*3+2]};
		real_hi vi[]={mVel[i*3], mVel[i*3+1], mVel[i*3+2]};
		real_lo ai[]={0,0,0};
		real_lo ji[]={0,0,0};

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
			real_hi r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
			real_hi rv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
			real_hi rinv = RSQRT(r2);
			rv *= 3./r2;
			rinv *= d_mass[j];
			real_hi rinv3 = rinv/r2;

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

//#pragma unroll 
	for(unsigned int i=1;i<nBodies;++i)
	{
		//float3 xi=mPos[i];
		real_hi xi[]={mPos[i*3], mPos[i*3+1], mPos[i*3+2]};
		real_hi vi[]={mVel[i*3], mVel[i*3+1], mVel[i*3+2]};
		real_lo ai[]={0,0,0};
		real_lo ji[]={0,0,0};

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
			real_hi r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
			//float r2 = dot(dx,dx);
			real_hi rv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
			//float rv = dot(dx,dv);
			real_hi rinv = RSQRT(r2);
			rv *= 3./r2;
			rinv *= d_mass[j];
			real_hi rinv3 = rinv/r2;

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
			real_hi r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
			//float r2 = dot(dx,dx);
			real_hi rv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
			//float rv = dot(dx,dv);
			real_hi rinv = RSQRT(r2);
			rv *= 3./r2;
			rinv *= d_mass[0];
			real_hi rinv3 = rinv/r2;

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

template<unsigned int pre, unsigned int nbod>
__global__ void gpu_hermite_integrator_kernel(double dT, double h)
{
	using namespace gpu_hermite_aux;

	ensemble &ens = gpu_hermite_ens;
	int sys = threadId();
	if(sys >= ens.nsys()) { return; }

	double    T = ens.time(sys);
	double Tend = T + dT;


	const unsigned int nData=3*nbod;
	typename pos_type<pre>::type mPos       [nData];
	typename pos_type<pre>::type mVel       [nData];
	typename acc_type<pre>::type mAcc       [nData];
	typename acc_type<pre>::type mJerk      [nData];
	typename pos_type<pre>::type mPosOld    [nData];
	typename pos_type<pre>::type mVelOld    [nData];
	typename acc_type<pre>::type mAccOld    [nData];
	typename acc_type<pre>::type mJerkOld   [nData];


	float s_mass[nbod];
	//const float s_mass[]={ens.mass(sys, 0), ens.mass(sys,1), ens.mass(sys,2), ens.mass(sys,3)};

	typename pos_type<pre>::type hh=h;
	typename acc_type<pre>::type dtby2=hh/2.;
	typename acc_type<pre>::type dtby3=hh/3.;
	typename acc_type<pre>::type dtby6=hh/6.;
	typename acc_type<pre>::type dt7by30=hh*7./30.;
	typename acc_type<pre>::type dtby7=hh*7.;

	//load data from global memory
	if(nbod>0)
	{
	mPos[0]=ens.x(sys,0);
	mPos[1]=ens.y(sys,0);
	mPos[2]=ens.z(sys,0);
	mVel[0]=ens.vx(sys,0);
	mVel[1]=ens.vy(sys,0);
	mVel[2]=ens.vz(sys,0);
	s_mass[0]=ens.mass(sys, 0);
	}
	if(nbod>1)
	{
	mPos[3]=ens.x(sys,1);
	mPos[4]=ens.y(sys,1);
	mPos[5]=ens.z(sys,1);
	mVel[3]=ens.vx(sys,1);
	mVel[4]=ens.vy(sys,1);
	mVel[5]=ens.vz(sys,1);
	s_mass[1]=ens.mass(sys, 1);
	}
	if(nbod>2)
	{
	mPos[6]=ens.x(sys,2);
	mPos[7]=ens.y(sys,2);
	mPos[8]=ens.z(sys,2);
	mVel[6]=ens.vx(sys,2);
	mVel[7]=ens.vy(sys,2);
	mVel[8]=ens.vz(sys,2);
	s_mass[2]=ens.mass(sys, 2);
	}
	if(nbod>3)
	{
	mPos[9]=ens.x(sys,3);
	mPos[10]=ens.y(sys,3);
	mPos[11]=ens.z(sys,3);
	mVel[9]=ens.vx(sys,3);
	mVel[10]=ens.vy(sys,3);
	mVel[11]=ens.vz(sys,3);
	s_mass[3]=ens.mass(sys, 3);
	}
	if(nbod>4)
	{
	mPos[12]=ens.x(sys,4);
	mPos[13]=ens.y(sys,4);
	mPos[14]=ens.z(sys,4);
	mVel[12]=ens.vx(sys,4);
	mVel[13]=ens.vy(sys,4);
	mVel[14]=ens.vz(sys,4);
	s_mass[4]=ens.mass(sys, 4);
	}
	if(nbod>5)
	{
	unsigned int idx = 15;
	for(unsigned int plid=5;plid<nbod;++plid)
		{	
		mPos[idx]=ens.x(sys,plid);
		mVel[idx]=ens.vx(sys,plid); ++idx;
		mPos[idx]=ens.y(sys,plid);
		mVel[idx]=ens.vy(sys,plid); ++idx;
		mPos[idx]=ens.z(sys,plid);
		mVel[idx]=ens.vz(sys,plid); ++idx;
		s_mass[plid]=ens.mass(sys, plid);
		}
	}

	if(pre==1)
	{
		if(nbod<5)
			UpdateAccJerk23<nbod>(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
		else
			UpdateAccJerkGeneral<nbod>(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
	}	
	else if(pre==2)
	{
// If not using shared memory as cache, no need bother	
//		float sPos       [nData];
//		float sVel       [nData];
//		copyArray<nData>(sPos,mPos);
//		copyArray<nData>(sVel,mVel);
		if(nbod<5)
			UpdateAccJerk23<nbod>(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
		else
			UpdateAccJerkGeneral<nbod>(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
	}
	else
	{
		float sPos       [nData];
		float sVel       [nData];
		copyArray<nData>(sPos,mPos);
		copyArray<nData>(sVel,mVel);
		if(nbod<5)
			UpdateAccJerk23<nbod>(&sPos[0], &sVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
		else
			UpdateAccJerkGeneral<nbod>(&sPos[0], &sVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
	}

	while(T<Tend)
	{

		copyArray<nData>(mPosOld,mPos);
		copyArray<nData>(mVelOld,mVel);
		copyArray<nData>(mAccOld,mAcc);
		copyArray<nData>(mJerkOld,mJerk);

	        if(T+h>Tend)
	          {	
			hh=Tend-T;
			dtby2=hh/2.;
			dtby3=hh/3.;
			dtby6=hh/6.;
			dt7by30=hh*7./30.;
			dtby7=hh*7.;
		  }
		predict<nData>(mPos,mVel,mAcc,mJerk, dtby2, dtby3, hh);

		if(pre==1)
		{
		if(nbod<5)
			UpdateAccJerk23<nbod>(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
		else
			UpdateAccJerkGeneral<nbod>(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
		}
		else if(pre==2)
		{
// If not using shared memory as cahce, don't bother
//		float sPos       [nData];
//		float sVel       [nData];
//		copyArray<nData>(sPos,mPos);
//		copyArray<nData>(sVel,mVel);
		if(nbod<5)
			UpdateAccJerk23<nbod>(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
		else
			UpdateAccJerkGeneral<nbod>(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
		}
		else
		{
		float sPos       [nData];
		float sVel       [nData];
		copyArray<nData>(sPos,mPos);
		copyArray<nData>(sVel,mVel);
		if(nbod<5)
			UpdateAccJerk23<nbod>(&sPos[0], &sVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
		else
			UpdateAccJerkGeneral<nbod>(&sPos[0], &sVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
		}

		correct<nData>(mPos,mVel,mAcc,mJerk, mPosOld,mVelOld,mAccOld,mJerkOld, dtby2, dtby6, dtby7, dt7by30);
		

		if(pre==1)
		{
		if(nbod<5)
			UpdateAccJerk23<nbod>(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
		else	
			UpdateAccJerkGeneral<nbod>(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
		} 
		else if(pre==2)
		{
// If not using shared memory as cahce, don't bother
//		float sPos       [nData];
//		float sVel       [nData];
//		copyArray<nData>(sPos,mPos);
//		copyArray<nData>(sVel,mVel);
		if(nbod<5)
			UpdateAccJerk23<nbod>(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
		else
			UpdateAccJerkGeneral<nbod>(&mPos[0], &mVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
		}
		else
		{
		float sPos       [nData];
		float sVel       [nData];
		copyArray<nData>(sPos,mPos);
		copyArray<nData>(sVel,mVel);
		if(nbod<5)
			UpdateAccJerk23<nbod>(&sPos[0], &sVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
		else
			UpdateAccJerkGeneral<nbod>(&sPos[0], &sVel[0], &mAcc[0], &mJerk[0], &s_mass[0]);
		}

		correct<nData>(mPos,mVel,mAcc,mJerk, mPosOld,mVelOld,mAccOld,mJerkOld, dtby2, dtby6, dtby7, dt7by30);

		T += hh;
		ens.nstep(sys)++;
	}
	ens.time(sys) = T;
	if(nbod>0)
	{
	ens.x(sys,0)=mPos[0];
	ens.y(sys,0)=mPos[1];
	ens.z(sys,0)=mPos[2];
	ens.vx(sys,0)=mVel[0];
	ens.vy(sys,0)=mVel[1];
	ens.vz(sys,0)=mVel[2];
	}
	if(nbod>1)
	{
	ens.x(sys,1)=mPos[3];
	ens.y(sys,1)=mPos[4];
	ens.z(sys,1)=mPos[5];
	ens.vx(sys,1)=mVel[3];
	ens.vy(sys,1)=mVel[4];
	ens.vz(sys,1)=mVel[5];
	}
	if(nbod>2)
	{
	ens.x(sys,2)=mPos[6];
	ens.y(sys,2)=mPos[7];
	ens.z(sys,2)=mPos[8];
	ens.vx(sys,2)=mVel[6];
	ens.vy(sys,2)=mVel[7];
	ens.vz(sys,2)=mVel[8];
	}
	if(nbod>3)
	{
	ens.x(sys,3)=mPos[9];
	ens.y(sys,3)=mPos[10];
	ens.z(sys,3)=mPos[11];
	ens.vx(sys,3)=mVel[9];
	ens.vy(sys,3)=mVel[10];
	ens.vz(sys,3)=mVel[11];
	}
	if(nbod>4)
	{
	ens.x(sys,4)=mPos[12];
	ens.y(sys,4)=mPos[13];
	ens.z(sys,4)=mPos[14];
	ens.vx(sys,4)=mVel[12];
	ens.vy(sys,4)=mVel[13];
	ens.vz(sys,4)=mVel[14];
	}
	if(nbod>5)
	{
	unsigned int idx = 15;
	for(unsigned int plid=5;plid<nbod;++plid)
		{	
		mPos[idx]=ens.x(sys,plid);
		mVel[idx]=ens.vx(sys,plid); ++idx;
		mPos[idx]=ens.y(sys,plid);
		mVel[idx]=ens.vy(sys,plid); ++idx;
		mPos[idx]=ens.z(sys,plid);
		mVel[idx]=ens.vz(sys,plid); ++idx;
		}
	}

}


template<>
void gpu_hermite_integrator<double,double>::integrate(gpu_ensemble &ens, double dT)
#include"hermite_gpu_integrator_body.cu"

template<>
void gpu_hermite_integrator<double,float>::integrate(gpu_ensemble &ens, double dT)
#include"hermite_gpu_integrator_body.cu"

template<>
void gpu_hermite_integrator<float,float>::integrate(gpu_ensemble &ens, double dT)
#include"hermite_gpu_integrator_body.cu"


} // end namespace hermite_gpu
} // end namespace swarm

