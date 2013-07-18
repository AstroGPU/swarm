#include <math.h>
#include <stdio.h>

#define max(a,b) ((a)<(b))?(b):(a)
#define sqr(x) (x)*(x)

typedef double REAL;
typedef int INT;

typedef void (*FCN)(INT,REAL,REAL*,REAL*,REAL*,INT*);

typedef void (*SOLFIX)(INT,REAL,REAL,REAL*, REAL*,INT,REAL*,INT*);

template<size_t nsd,size_t nmd>
void coef(INT ns, REAL* C, REAL* B, REAL* BC, REAL (&AA)[nsd][nsd], REAL (&E)[nsd][nsd+nmd], REAL* SM, REAL* AM, REAL hStep);

template<size_t ndgl,size_t nsd, size_t nmd>
void startb(FCN fcn, INT N, REAL X, REAL h,REAL* P, REAL* Q, INT ns, REAL* FS, REAL* PS, REAL (&ZQ)[ndgl][nsd], REAL (&E)[nsd][nsd+nmd], REAL* YH, REAL* SM, REAL* AM, REAL* F, REAL* C, REAL* RPAR, INT* IPAR);

template<size_t ndgl, size_t nsd>
void rknife(FCN fcn, INT N, INT ns, REAL X, REAL* Q, REAL* P,REAL (&AA)[nsd][nsd], REAL* C, REAL* QQ, REAL (&ZQ)[ndgl][nsd],REAL* F, REAL& dyno, REAL* RPAR, INT* IPAR);
/**
* C  SOLVES SECOND ORDER ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM
C                       Q'' = F(X,Q)
C  BASED ON THE SYMPLECTIC AND SYMMETRIC GAUSS (IRK) METHODS
C  DESCRIBED IN SECTIONS II.1, VIII.6 OF THE BOOK:
C
C      E. HAIRER, C. LUBICH, G. WANNER, GEOMETRIC NUMERICAL INTEGRATION,
C         STRUCTURE-PRESERVING ALGORITHMS FOR ODES.
C         SPRINGER SERIES IN COMPUT. MATH. 31, SPRINGER 2002.
C
C  AND IN THE PUBLICATION
C
C      E. HAIRER, M. HAIRER, GNI-CODES - MATLAB PROGRAMS FOR
C         GEOMETRIC NUMERICAL INTEGRATION.
C
C  INPUT..
C     N           DIMENSION OF Q AND F(X,Q) 
C
C     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING F(X,Q):
C                    SUBROUTINE FCN(N,X,Q,F,RPAR,IPAR)
C                    REAL*8 Q(N),F(N)
C                    F(1)=...   ETC.
C
C     NSTEP       NUMBER OF INTEGRATION STEPS
C                    CONSTANT STEP SIZE, H=(XEND-X)/NSTEP
C
C     X           INITIAL X-VALUE
C     P(N)        INITIAL VELOCITY VECTOR
C     Q(N)        INITIAL POSITION VECTOR
C     XEND        FINAL X-VALUE
C
C     METH        NUMBER OF STAGES OF THE GAUSS METHOD
C                    FOR THE MOMENT ONLY POSSIBLE VALUES: 2,4,6.
C
C     SOLFIX      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
C                 NUMERICAL SOLUTION DURING INTEGRATION. 
C                 IF IOUT=1, IT IS CALLED AFTER EVERY STEP.
C                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. 
C                    SUBROUTINE SOLFIX (NR,XOLD,X,P,Q,N,IRTRN,RPAR,IPAR)
C                    DOUBLE PRECISION X,Y(N),CONT(LRC)
C                      ....  
C                 SOLFIX FURNISHES THE SOLUTION "Q,P" AT THE NR-TH
C                    GRID-POINT "X" (INITIAL VALUE FOR NR=0).
C                 "XOLD" IS THE PRECEEDING GRID-POINT.
C                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
C                    IS SET <0, RETURN TO THE CALLING PROGRAM.
C     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLFIX:
C                    IOUT=0: SUBROUTINE IS NEVER CALLED
C                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT.
C
C     RPAR(LR)    REAL PARAMETER ARRAY; LR MUST BE AT LEAST LR=10
C                    RPAR(1),...,RPAR(10) SERVE AS PARAMETERS FOR
C                    THE CODE. FURTHER VALUES CAN BE USED FOR DEFINING
C                    PARAMETERS IN THE PROBLEM
C     IPAR(LI)    INTEGER PARAMETER ARRAY; LI MUST BE AT LEAST LI=10
C                    IPAR(1),...,IPAR(10) SERVE AS PARAMETERS FOR
C                    THE CODE. FURTHER VALUES CAN BE USED FOR DEFINING
C                    PARAMETERS IN THE PROBLEM
C
C  OUTPUT..
C     P(N)        SOLUTION (VELOCITY) AT XEND
C     Q(N)        SOLUTION (POSITION) AT XEND
C-----------------------------------------------------------------------
C     SOPHISTICATED SETTING OF PARAMETERS 
C-----------------------------------------------------------------------
C    RPAR(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
C    IPAR(1)   NITMAX, MAXIMAL NUMER OF FIXED POINT ITERAT., DEFAULT 50
C-----------------------------------------------------------------------
*/
inline void gni_irk2(INT N, FCN fcn,INT nStep,REAL X,REAL* P, REAL* Q, REAL XEnd,INT method, SOLFIX solfix,bool IOUT, REAL* RPAR, INT* IPAR)
{
	const static INT nsd = 6, nmd = 3, ndgl = 2; // ndgl = N
	REAL F[ndgl*nsd], YH[ndgl], QQ[ndgl];
	REAL C[nsd],AA[nsd][nsd],E[nsd][nsd+nmd],B[nsd],BC[nsd];
	REAL SM[nmd],AM[nsd+nmd];
	REAL FS[ndgl],PS[ndgl], ZQ[ndgl][nsd];

	REAL uround;
	if (RPAR[0] == 0.0)
		uround = 1e-16;
	else
		uround = RPAR[0];
	int nitmax;
	if (IPAR[0] == 0)
		nitmax = 50;
	else
		nitmax = IPAR[0];
	
	INT ns = method;
	REAL h = (XEnd - X)/nStep;
	
	coef<nsd,nmd>(ns,C,B,BC,AA,E,SM,AM,h);

	if (IOUT) solfix(0,X,X,P,Q,N,RPAR,IPAR);

	fcn(N,X,Q,FS,RPAR,IPAR);
	for(int is = 0; is<ns; is++)
	{
		REAL FAC = C[is]*C[is]/2.0;
		for (int i = 0; i<N; i++)
			ZQ[i][is] = C[is]*P[i] + FAC*FS[i];

	}
	for(int i = 0; i<N; i++)
	{
		PS[i] = P[i];
		
	}
	// main loop
	REAL dynold, dyno;
	for(int istep = 0; istep < nStep; istep++)
	{
		//printf("%d..",istep);
		if (istep > 0) 
			startb<ndgl,nsd,nmd>(fcn,N,X,h,P,Q,ns,FS,PS,ZQ,E,YH,SM,AM,F,C,RPAR,IPAR);
		// fixed point iteration
		INT niter = 0;
		dynold = 0.0;
		dyno = 1.0;
		while (dyno > uround)
		{
			rknife<ndgl,nsd>(fcn,N,ns,X,Q,P,AA,C,QQ,ZQ,F,dyno,RPAR,IPAR);
			
			niter++;
			if ((dynold < dyno) && (dyno < 10*uround)) 
				break;
			if (niter >= nitmax)
			{
				printf("no convergence of iteration: %f\n", dyno);
				return;
			}
			dynold = dyno;
		}
		// update of the solution
		X = X + h;
		for(int i = 0; i < N; i++)
		{
			REAL sum = 0.0;
			for(int is = 0; is<ns; is++)
				sum += F[i+is*N]*BC[is];
			Q[i] += h*P[i] + sum;

			sum = 0.0;
			for(int is = 0; is<ns; is++)
				sum += F[i+is*N]*B[is];
			P[i] += sum;
			
			//printf("q%d = %f, p%d = %f ", i, Q[i], i, P[i]);
		}
		//printf("\n");
		if (IOUT) solfix(istep,X-h,X,P,Q,N,RPAR,IPAR);
	}
}
template<size_t ndgl,size_t nsd, size_t nmd>
void startb(FCN fcn, INT N, REAL X, REAL h,REAL* P, REAL* Q, INT ns, REAL* FS, REAL* PS, REAL (&ZQ)[ndgl][nsd], REAL (&E)[nsd][nsd+nmd], REAL* YH, REAL* SM, REAL* AM, REAL* F, REAL* C, REAL* RPAR, INT* IPAR)
{
	INT ns1 = ns;
	INT ns2 = ns + 1;
	INT nsm = ns + nmd - 1;
	for( int i = 0; i < N; i++)
	{
		REAL sav = 0.0;
		for(int js = 0; js<ns; js++)
			sav += AM[js]*ZQ[i][js];
		YH[i] = sav + AM[ns1]*PS[i] + AM[ns2]*P[i]+Q[i];
		for(int is=0; is < ns; is++)
		{
			REAL sav = 0.0;
			for(int js = 0; js < ns; js++)
				sav += E[is][js]*F[i+js*N];
			ZQ[i][is] = sav + E[is][ns1]*FS[i];
		}
	}
	fcn(N,X+h,Q,FS,RPAR,IPAR);
	fcn(N,X+h*SM[nmd-1],YH,F,RPAR,IPAR);
	for(int i = 0; i<N; i++)
	{
		PS[i] = P[i];
		for (int is = 0; is < ns;is++)
		{
			ZQ[i][is] += E[is][ns2]*FS[i]+E[is][nsm]*F[i] + C[is]*P[i]; 
		}
	}
}
template<size_t ndgl, size_t nsd>
void rknife(FCN fcn, INT N, INT ns, REAL X, REAL* Q, REAL* P,REAL (&AA)[nsd][nsd], REAL* C, REAL* QQ, REAL (&ZQ)[ndgl][nsd],REAL* F, REAL& dyno, REAL* RPAR, INT* IPAR)
{
	for(int js=0; js<ns; js++)
	{
		for(int j=0; j<N; j++)
			QQ[j] = Q[j] + ZQ[j][js];
		fcn(N,X+C[js],QQ,F+js*N,RPAR,IPAR); 
	}
	dyno = 0.0;
	for(int i = 0; i<N; i++)
	{
		REAL dnom = max(1e-1,abs(Q[i]));
		for(int is = 0; is<ns; is++)
		{
			REAL sum = C[is]*P[i];
			for(int js=0; js<ns; js++)
				sum += AA[is][js]*F[i + js*N];
			dyno += sqr((sum - ZQ[i][is])/dnom);
			ZQ[i][is] = sum;
		}

	}
	dyno = sqrt(dyno/(ns*N));

}
template<size_t nsd,size_t nmd>
void coef(INT ns, REAL* C, REAL* B, REAL* BC, REAL (&AA)[nsd][nsd], REAL (&E)[nsd][nsd+nmd], REAL* SM, REAL* AM, REAL hStep)
{
	
	if (ns == 2)
	{
		 C[0]= 0.21132486540518711775;
         C[1]= 0.78867513459481288225;
         B[0]= 0.50000000000000000000;
         B[1]= 0.50000000000000000000;
         BC[0]= 0.39433756729740644113;
         BC[1]= 0.10566243270259355887;
         AA[0][0]= 0.41666666666666666667e-1;
         AA[0][1]=-0.19337567297406441127e-1;
         AA[1][0]= 0.26933756729740644113e+0;
         AA[1][1]= 0.41666666666666666667e-1;
         E[0][0]=-0.28457905077110526160e-02;
         E[0][1]=-0.63850024471784160410e-01;
         E[0][2]= 0.48526095198694517563e-02;
         E[0][3]= 0.11305688530429939012e+00;
         E[0][4]=-0.28884580475413403312e-01;
         E[1][0]= 0.41122751744511433137e-01;
         E[1][1]=-0.18654814888622834132e+00;
         E[1][2]=-0.18110185277445209332e-01;
         E[1][3]= 0.36674109449368040786e+00;
         E[1][4]= 0.10779872188955481745e+00;
         SM[0]= 0.00000000000000000000e+00;
         SM[1]= 0.10000000000000000000e+01;
         SM[2]= 0.16000000000000000000e+01;
         AM[0]= 0.25279583039343438291e+02;
         AM[1]=-0.86907830393434382912e+01;
         AM[2]=-0.80640000000000000000e+00;
         AM[3]= 0.29184000000000000000e+01;
         AM[4]= 0.00000000000000000000e+00;
	}
	if (ns==4)
	{
         C[0]= 0.69431844202973712388e-01;
         C[1]= 0.33000947820757186760e+00;
         C[2]= 0.66999052179242813240e+00;
         C[3]= 0.93056815579702628761e+00;
         B[0]= 0.17392742256872692869e+00;
         B[1]= 0.32607257743127307131e+00;
         B[2]= 0.32607257743127307131e+00;
         B[3]= 0.17392742256872692869e+00;
         BC[0]= 0.16185132086231030665e+00;
         BC[1]= 0.21846553629538057030e+00;
         BC[2]= 0.10760704113589250101e+00;
         BC[3]= 0.12076101706416622036e-01;
         AA[0][0]= 0.40381914508467311298e-02;
         AA[0][1]=-0.32958609449446961650e-02;
         AA[0][2]= 0.26447829520668538006e-02;
         AA[0][3]=-0.97672296325588161023e-03;
         AA[1][0]= 0.43563580902396261254e-01;
         AA[1][1]= 0.13818951406296126013e-01;
         AA[1][2]=-0.43401341944349953440e-02;
         AA[1][3]= 0.14107297391595337720e-02;
         AA[2][0]= 0.10586435263357640763e+00;
         AA[2][1]= 0.10651836096505307395e+00;
         AA[2][2]= 0.13818951406296126013e-01;
         AA[2][3]=-0.17580153590805494993e-02;
         AA[3][0]= 0.14879849619263780300e+00;
         AA[3][1]= 0.19847049885237718995e+00;
         AA[3][2]= 0.81671359795877570687e-01;
         AA[3][3]= 0.40381914508467311298e-02;
         E[0][0]=-0.21272768296134340207e-1;
         E[0][1]= 0.11059138674756969912e-1;
         E[0][2]= 0.38999255049973564023e-2;
         E[0][3]=-0.43986226789008967612e-1;
         E[0][4]= 0.13581590305438849621e-1;
         E[0][5]= 0.39922421675314269059e-1;
         E[0][6]=-0.79369058065113002021e-3;
         E[1][0]=-0.75671119283734809953e-02;
         E[1][1]= 0.10209394000843457002e-01;
         E[1][2]=-0.12880197817980892596e-01;
         E[1][3]=-0.56381316813776501277e-01;
         E[1][4]= 0.37440782682669799960e-02;
         E[1][5]= 0.11522469441011273193e+00;
         E[1][6]= 0.21035877343246316334e-02;
         E[2][0]=-0.39890571772473709759e+00;
         E[2][1]= 0.26819725655216894347e+00;
         E[2][2]=-0.82551711648854471247e-01;
         E[2][3]=-0.85516559106259630212e+00;
         E[2][4]= 0.24433810515772642570e+00;
         E[2][5]= 0.10234155624049009806e+01;
         E[2][6]= 0.25115745967236579242e-01;
         E[3][0]=-0.40964796048052939224e+00;
         E[3][1]= 0.29949323098224574487e+00;
         E[3][2]=-0.13867460566101912494e+00;
         E[3][3]=-0.98859300714628940382e+00;
         E[3][4]= 0.24671351779481625627e+00;
         E[3][5]= 0.12912760231350872304e+01;
         E[3][6]= 0.13241134766742798418e+00;
         SM[0]= 0.00000000000000000000e+00;
         SM[1]= 0.10000000000000000000e+01;
         SM[2]= 0.16500000000000000000e+01;
         AM[0]= 0.10806374869244001787e+04;
         AM[1]=-0.66008818661284690206e+03;
         AM[2]= 0.61810154357557529566e+03;
         AM[3]=-0.31341427826212857229e+03;
         AM[4]=-0.10187174765625000000e+02;
         AM[5]= 0.31173050390625000000e+02;
         AM[6]= 0.00000000000000000000e+00;
	}
	REAL hstep2 = hStep*hStep;
	for (int i = 0; i<ns; i++)
	{
		B[i] *= hStep;
		BC[i] *= hstep2;
		C[i] *= hStep;
		for(int j = 0; j<ns; j++)
		{
			AA[i][j] *=hstep2;
			E[i][j] *=hstep2;
		}
		
	}
	for(int i1 = 0; i1<nmd; i1++)
	{
		for(int i2 = 0; i2<ns; i2++)
		{
			E[i2][i1+ns] *=hstep2;
		}
		AM[ns+i1] *=hStep;
	}

}


