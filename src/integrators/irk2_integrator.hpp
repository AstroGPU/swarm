/*************************************************************************
 * Copyright (C) 2013 by Thien Nguyen and the Swarm-NG Development Team  *
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

/*! \file irk2_integrator.hpp
 *   \brief Defines and implements \ref swarm::gpu::bppt::irk2_integrator class - the 
 *          GPU implementation of Implicit Runge Kutta described in 
 *          E. HAIRER, M. HAIRER, GNI-CODES - MATLAB PROGRAMS FOR
 *          GEOMETRIC NUMERICAL INTEGRATION..
 *
 */

#include "swarm/common.hpp"
#include "swarm/gpu/bppt.hpp"

namespace swarm { namespace gpu { namespace bppt {

/*! GPU implementation of Implicit Runge Kutta
 * \ingroup integrators
 *
 */

template< class Monitor , template<class T> class Gravitation > //Gravitation should be GravitationAcc defined in gravitation_acc.hpp
class irk2: public integrator {
	typedef integrator base;
	typedef Monitor monitor_t;
	typedef typename monitor_t::params mon_params_t;
	
private:
	double _time_step;
	
	int ns; // method
	
	mon_params_t _mon_params;

public: //! Construct for class hermite integrator
	irk2(const config& cfg): base(cfg),_time_step(0.001), _mon_params(cfg) {
		_time_step =  cfg.require("time_step", 0.0);
		ns = cfg.optional("method",4);
		
	}

	virtual void launch_integrator() {
		launch_templatized_integrator(this);
	}
	
	
        //! Define the number of thread per system
	template<class T>
	static GENERIC int thread_per_system(T compile_time_param){
		const int grav = Gravitation<T>::thread_per_system();
		const int moni = Monitor::thread_per_system(compile_time_param);
		return max( grav, moni);
	}

        //! Define the amount of shared memory per system
	template<class T>
	static GENERIC int shmem_per_system(T compile_time_param){
		const int grav = Gravitation<T>::shmem_per_system();
		const int moni = Monitor::shmem_per_system(compile_time_param);
		return max( grav, moni);
	}	

        //! Convert internal coordinates to std coordinates
        GPUAPI void convert_internal_to_std_coord() {} 
        //! Convert std coordinates to internal coordinates
        GPUAPI void convert_std_to_internal_coord() {}  

	template<class T>
	__device__ void kernel(T compile_time_param){

		if(sysid()>=_dens.nsys()) return;
		
		// References to Ensemble and Shared Memory
		typedef Gravitation<T> Grav;
		ensemble::SystemRef sys = _dens[sysid()];
		
		typedef typename Monitor::shared_data data_t;
		monitor_t montest(_mon_params,sys,*_log, *((data_t*) system_shared_data_pointer(this,compile_time_param))) ;
				
		typedef typename Grav::shared_data grav_t;
		Grav calcForces(sys,*( (grav_t*) system_shared_data_pointer(this,compile_time_param) ) );
		
		

		// Local variables
		const int nbod = T::n;
		// Body number
		const int b = thread_body_idx(nbod);
		// Component number
		const int c = thread_component_idx(nbod);
		
		const int nsd = 6, nmd = 3;
		// the following coefficients should be declared as __constant__, but I dont know how 
		double C[nsd],AA[nsd][nsd],E[nsd][nsd+nmd],B[nsd],BC[nsd],SM[nmd],AM[nsd+nmd];
		
		//if (thread_in_system() == 0)
		coef<nsd,nmd>(ns,C,B,BC,AA,E,SM,AM,_time_step);
		//__syncthreads();
		
		double F[nsd], YH, QQ, FS, PS, ZQ[nsd];
		
		const double uround = 1e-16;
		
		// local variables
		montest.init( thread_in_system() );
		__syncthreads();

		// local information per component per body
		double pos = 0.0, vel = 0.0 , acc0 = 0.0 , c_pos, c_vel; // P = vel; Q = pos;
		if( (b < nbod) && (c < 3) )
			{ pos = sys[b][c].pos(); vel = sys[b][c].vel(); }


		////////// INTEGRATION //////////////////////

		/// Calculate acceleration
		calcForces(thread_in_system(),b,c,pos,vel,acc0);
		
		if( (b < nbod) && (c < 3) )
		{
			FS = acc0; 
			for (int is = 0; is<ns; is++)
				ZQ[is] = C[is]*sys[b][c].vel() + 0.5*C[is]*C[is]*FS;
				
			PS = vel;
		}
		
		//typedef DoubleCoalescedStruct<SHMEM_CHUNK_SIZE> shared_para_t[3]; // shared data for nit, dynold, dyno in nonlinear solver
		//shared_para_t& shared_para = * (shared_para_t*) system_shared_data_pointer(this,compile_time_param) ;
		
                for(int iter = 0 ; (iter < _max_iterations) && sys.is_active() ; iter ++ ) 
		{
			double h = _time_step;
                        
			if( sys.time() + h > _destination_time ) {
				h = _destination_time - sys.time();
			}
			c_pos = sys[b][c].pos();
                        c_vel = sys[b][c].vel();
			// startb
			if (iter > 0)
			{
				
				double sav = 0.0;
				for(int js = 0; js<ns; js++)
					sav += AM[js]*ZQ[js];
				YH = sav + AM[ns]*PS + AM[ns+1]*vel+pos;
				for(int is=0; is < ns; is++)
				{
					sav = 0.0;
					for(int js = 0; js < ns; js++)
						sav += E[is][js]*F[js];
					ZQ[is] = sav + E[is][ns]*FS;
				}
				calcForces(thread_in_system(),b,c,pos,vel,FS);
				calcForces(thread_in_system(),b,c,YH,vel,F[0]);
				PS = vel;
				for (int is = 0; is < ns;is++)
				{
					ZQ[is] += E[is][ns+1]*FS+E[is][ns + nmd - 1]*F[0] + C[is]*vel; 
				}
				
			}// end of startb			
			
			/// fixed point iteration
			if (thread_in_system() == 0)
			{
				sys.attribute(0) = 0;//int nit = 0;
				sys.attribute(1) = 0.0;//double dynold = 0.0;
				sys.attribute(2) = 1.0;//double dyno = 1.0;
			}
			__syncthreads();	
			while (sys.attribute(2) > uround)
			{
				/// rknife
				if (thread_in_system() == 0)
				{
					sys.attribute(2) = 0.0;
				}
				__syncthreads();
				
				for(int js=0; js<ns; js++)
				{
					QQ = pos + ZQ[js];
					calcForces(thread_in_system(),b,c,QQ,vel,F[js]); 
				}
				
				double dnom = max(1e-1,abs(pos));
                                double ss = 0;
				for(int is = 0; is<ns; is++)
				{
					double sum = C[is]*vel;
					for(int js=0; js<ns; js++)
						sum += AA[is][js]*F[js];
					ss += (sum - ZQ[is])*(sum - ZQ[is])/(dnom*dnom);
                                       
					ZQ[is] = sum;
				}
				atomicAdd(&sys.attribute(2),ss);
                                
				__syncthreads();

				if (thread_in_system() == 0)
				{
					sys.attribute(2) = sqrt(sys.attribute(2)/(ns*3*nbod));
					// end of rknife
                                        
					sys.attribute(0) = sys.attribute(0) + 1;
                                                                                
                                        //lprintf(*_log,"%d %d %d: %f %f\n", iter,b,c, sys.attribute(0),sys.attribute(2));
                                        
                                        if (sys.attribute(0) >= 50)
                                        {
                                                lprintf(*_log,"no convergence of iteration: %f\n", sys.attribute(0));
                                                sys.set_inactive();
                                                
                                        }
                                                                  
				}
				__syncthreads();
                                
                                if (sys.attribute(0) >= 50) break;
                                                                
                                if ((sys.attribute(1) < sys.attribute(2)) && (sys.attribute(2) < 10*uround)) 
                                   break;
                                
                                __syncthreads();
                                
                                if (thread_in_system() == 0)
                                {
                                  sys.attribute(1) = sys.attribute(2);
                                }
                                
                                
			}// end of fixed point iteration.
			
			// update solution
			
			double sum = 0.0;
			for(int is = 0; is<ns; is++)
				sum += F[is]*BC[is];
			pos = pos + h*vel + sum;

			sum = 0.0;
			for(int is = 0; is<ns; is++)
				sum += F[is]*B[is];
			vel = vel + sum;
			
                        //lprintf(*_log,"%d %d %d: %f %f\n", iter, b,c, pos,vel);
                        /// Finalize the step
			if( (b < nbod) && (c < 3) )
				{ sys[b][c].pos() = pos; sys[b][c].vel() = vel; }
			if( thread_in_system()==0 ) 
				sys.time() += h;
			
			/// Monitor collision-detection
//                         montest.storeCurrentStat(b,c,c_pos,c_vel);
//                         montest(thread_in_system(), b, c);  
//                         __syncthreads();
			
			if( sys.is_active() && thread_in_system()==0 )  {
			    if( sys.time() >= _destination_time ) 
			    {	sys.set_inactive(); }
			}

			__syncthreads();


		}

	}
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
template<size_t nsd,size_t nmd>
GPUAPI void coef(int ns, double* C, double* B, double* BC, double (&AA)[nsd][nsd], double (&E)[nsd][nsd+nmd], double* SM, double* AM, double hStep)
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
	double hstep2 = hStep*hStep;
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


};

} } } // end namespace bppt :: integrators :: swarm
