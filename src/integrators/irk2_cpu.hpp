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

/*! \file irk2_cpu.hpp
 *   \brief Defines and implements \ref swarm::cpu::irk2_cpu class - the 
 *          CPU implementation of PEC2 Hermite integrator.
 *
 */

#ifdef _OPENMP
#include <omp.h>
#endif


#include "swarm/common.hpp"
#include "swarm/integrator.hpp"
#include "swarm/plugin.hpp"

namespace swarm { namespace cpu {
/*! CPU implementation of implicit Runge-Kutta integrator
 *
 * \ingroup integrators
 *
 *   This is used as a reference implementation to
 *   test other GPU implementations of other integrators
 *   
 *   This integrator can be used as an example of CPU integrator
 *
 */
template< class Monitor >
class irk2_cpu : public integrator {
        typedef integrator base;
        typedef Monitor monitor_t;
        typedef typename monitor_t::params mon_params_t;
        private:
        double _time_step;
        mon_params_t _mon_params;

public:  //! Construct for hermite_cpu class
        irk2_cpu(const config& cfg): base(cfg),_time_step(0.001), _mon_params(cfg) {
                _time_step =  cfg.require("time_step", 0.0);
        }

        virtual void launch_integrator() {
                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
                for(int i = 0; i < _ens.nsys(); i++){
                        integrate_system(_ens[i]);
                }
        }

        //! defines inner product of two arrays
        inline static double inner_product(const double a[3],const double b[3]){
                return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
        }

        //! Calculate the force field. 
        void calcForces(ensemble::SystemRef& sys,double* pos, double* acc){
               
          const int nbod = sys.nbod();
               
               for(int b = 0; b < 3*nbod; b++)   
                  acc[b] = 0;

                /// Loop through all pairs
                for(int i=0; i < nbod-1; i++) 
                    for(int j = i+1; j < nbod; j++) 
                    {

                        double dx[3] = { pos[3*j]- pos[3*i],
                                pos[3*j+1]- pos[3*i+1],
                                pos[3*j+2]- pos[3*i+2],
                        };
                        
                        /// Calculated the magnitude
                        double r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2] * dx[2];
                        double rinv = 1 / ( sqrt(r2) * r2 ) ;
                        
                        /// Update acc/jerk for i
                        const double scalar_i = +rinv*sys[j].mass();
                        for(int c = 0; c < 3; c++) {
                                acc[3*i+c] += dx[c]* scalar_i;
                                
                        }

                        /// Update acc/jerk for j
                        const double scalar_j = -rinv*sys[i].mass();
                        for(int c = 0; c < 3; c++) {
                                acc[3*j+c] += dx[c]* scalar_j;
                                
                        }
                }
        }

        //! Integrate ensembles
        void integrate_system(ensemble::SystemRef sys){
                const int nbod = sys.nbod();
                const static int nsd = 6, nmd = 3, ndgl = 3*nbod;
                double F[ndgl*nsd], YH[ndgl], QQ[ndgl];
                double C[nsd],AA[nsd][nsd],E[nsd][nsd+nmd],B[nsd],BC[nsd];
                double SM[nmd],AM[nsd+nmd];
                double FS[ndgl],PS[ndgl], ZQ[ndgl][nsd];

                double uround=1e-16;
                
                int nitmax = 50;
                int ns = 6;
                double h = _time_step;
                int N = ndgl;
                
                coef<nsd,nmd>(ns,C,B,BC,AA,E,SM,AM,h);
                
                double* Q = new double[N];
                double* P = new double[N];
                for (int i = 0; i < nbod; i++)
                {
                  for (int j = 0; j < 3; j++)
                  {
                    Q[3*i+j] = sys[i][j].pos();
                    P[3*i+j] = sys[i][j].vel();
                  }
                }
                calcForces(sys,Q,FS);
                
                for(int is = 0; is<ns; is++)
                {
                        double FAC = C[is]*C[is]/2.0;
                        for (int i = 0; i<N; i++)
                                ZQ[i][is] = C[is]*P[i] + FAC*FS[i];

                }
                for(int i = 0; i<N; i++)
                {
                        PS[i] = P[i];
                        
                }
                // main loop
                double dynold, dyno;

                monitor_t montest (_mon_params,sys,*_log);


                for(int iter = 0 ; (iter < _max_iterations) && sys.is_active() ; iter ++ ) {
                
                        if( sys.time() + h > _destination_time ) {
                                h = _destination_time - sys.time();
                                coef<nsd,nmd>(ns,C,B,BC,AA,E,SM,AM,h);
                        }
                        
                        
                        if (iter > 0) 
                        {
                          int ns1 = ns;
                          int ns2 = ns + 1;
                          int nsm = ns + nmd - 1;
                          for( int i = 0; i < N; i++)
                          {
                                  double sav = 0.0;
                                  for(int js = 0; js<ns; js++)
                                          sav += AM[js]*ZQ[i][js];
                                  YH[i] = sav + AM[ns1]*PS[i] + AM[ns2]*P[i]+Q[i];
                                  for(int is=0; is < ns; is++)
                                  {
                                          double sav = 0.0;
                                          for(int js = 0; js < ns; js++)
                                                  sav += E[is][js]*F[i+js*N];
                                          ZQ[i][is] = sav + E[is][ns1]*FS[i];
                                  }
                          }
                          calcForces(sys,Q,FS);
                          calcForces(sys,YH,F);
                          for(int i = 0; i<N; i++)
                          {
                                  PS[i] = P[i];
                                  for (int is = 0; is < ns;is++)
                                  {
                                          ZQ[i][is] += E[is][ns2]*FS[i]+E[is][nsm]*F[i] + C[is]*P[i]; 
                                  }
                          }
                          
                        }
                        // fixed point iteration
                        int niter = 0;
                        dynold = 0.0;
                        dyno = 1.0;
                        while (dyno > uround)
                        {
                                for(int js=0; js<ns; js++)
                                {
                                        for(int j=0; j<N; j++)
                                                QQ[j] = Q[j] + ZQ[j][js];
                                        calcForces(sys,QQ,F+js*N); 
                                }
                                dyno = 0.0;
                                for(int i = 0; i<N; i++)
                                {
                                        double dnom = std::max(1e-1,std::abs(Q[i]));
                                        for(int is = 0; is<ns; is++)
                                        {
                                                double sum = C[is]*P[i];
                                                for(int js=0; js<ns; js++)
                                                        sum += AA[is][js]*F[i + js*N];
                                                dyno += (sum - ZQ[i][is])*(sum - ZQ[i][is])/(dnom*dnom);
                                                ZQ[i][is] = sum;
                                        }

                                }
                                dyno = sqrt(dyno/(ns*N));
                                
                                
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
                        
                        for(int b = 0; b < nbod; b++)
                        for(int c = 0; c < 3; c++)
                        {
                                int i = 3*b + c;
                                double sum = 0.0;
                                for(int is = 0; is<ns; is++)
                                        sum += F[i+is*N]*BC[is];
                                Q[i] += h*P[i] + sum;
                                
                                sys[b][c].pos() = Q[i];

                                sum = 0.0;
                                for(int is = 0; is<ns; is++)
                                        sum += F[i+is*N]*B[is];
                                P[i] += sum;
                                
                                sys[b][c].vel() = P[i];
                                                   
                        }
                                                             
                        sys.time() += h;
                        
                        if( sys.is_active() )  {
                                montest(0);
                                if( sys.time() > _destination_time - 1e-12) 
                                        sys.set_inactive();
                        }

                }
        }
template<size_t nsd,size_t nmd>
void coef(int ns, double* C, double* B, double* BC, double (&AA)[nsd][nsd], double (&E)[nsd][nsd+nmd], double* SM, double* AM, double hStep)
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
        if (ns == 6)
        {
          C[0]= 0.33765242898423986094e-01;
          C[1]= 0.16939530676686774317e+00;
          C[2]= 0.38069040695840154568e+00;
          C[3]= 0.61930959304159845432e+00;
          C[4]= 0.83060469323313225683e+00;
          C[5]= 0.96623475710157601391e+00;
          B[0]= 0.85662246189585172520e-01;
          B[1]= 0.18038078652406930378e+00;
          B[2]= 0.23395696728634552369e+00;
          B[3]= 0.23395696728634552369e+00;
          B[4]= 0.18038078652406930378e+00;
          B[5]= 0.85662246189585172520e-01;
          BC[0]= 0.82769839639769234611e-01;
          BC[1]= 0.14982512785597570103e+00;
          BC[2]= 0.14489179419935320895e+00;
          BC[3]= 0.89065173086992314743e-01;
          BC[4]= 0.30555658668093602753e-01;
          BC[5]= 0.28924065498159379092e-02;
          AA[0][0]= 0.90625420195651151857e-03;
          AA[0][1]=-0.72859711612531400024e-03;
          AA[0][2]= 0.79102695861167691135e-03;
          AA[0][3]=-0.70675390218535384182e-03;
          AA[0][4]= 0.45647714224056921122e-03;
          AA[0][5]=-0.14836147050330408643e-03;
          AA[1][0]= 0.11272367531794365387e-01;
          AA[1][1]= 0.39083482447840698486e-02;
          AA[1][2]=-0.14724868010943911900e-02;
          AA[1][3]= 0.10992669056588431310e-02;
          AA[1][4]=-0.67689040729401428165e-03;
          AA[1][5]= 0.21677950347174141516e-03;
          AA[2][0]= 0.30008019623627547434e-01;
          AA[2][1]= 0.36978289259468146662e-01;
          AA[2][2]= 0.65490339168957822692e-02;
          AA[2][3]=-0.16615098173008262274e-02;
          AA[2][4]= 0.84753461862041607649e-03;
          AA[2][5]=-0.25877462623437421721e-03;
          AA[3][0]= 0.49900269650650898941e-01;
          AA[3][1]= 0.82003427445271620462e-01;
          AA[3][2]= 0.54165111295060067982e-01;
          AA[3][3]= 0.65490339168957822692e-02;
          AA[3][4]=-0.11352871017627472322e-02;
          AA[3][5]= 0.28963081055952389031e-03;
          AA[4][0]= 0.68475836671617248304e-01;
          AA[4][1]= 0.11859257878058808400e+00;
          AA[4][2]= 0.10635984886129551097e+00;
          AA[4][3]= 0.47961474042181382443e-01;
          AA[4][4]= 0.39083482447840698486e-02;
          AA[4][5]=-0.34600839001342442657e-03;
          AA[5][0]= 0.79729071619449992615e-01;
          AA[5][1]= 0.14419100392702230613e+00;
          AA[5][2]= 0.13628542646896576408e+00;
          AA[5][3]= 0.81956586217401900627e-01;
          AA[5][4]= 0.23736460480774324642e-01;
          AA[5][5]= 0.90625420195651151857e-03;
          E[0][0]=-0.16761132335280609813e-01;
          E[0][1]= 0.10201050166615899799e-01;
          E[0][2]=-0.58593121685075943100e-02;
          E[0][3]=-0.11907383391366998251e-03;
          E[0][4]= 0.10615611118132982241e-01;
          E[0][5]=-0.30692054230989138447e-01;
          E[0][6]= 0.10615182045216224925e-01;
          E[0][7]= 0.22586707045496892369e-01;
          E[0][8]=-0.16931992776201068110e-04;
          E[1][0]= 0.10671755276327262128e-01;
          E[1][1]=-0.51098203653251450913e-02;
          E[1][2]= 0.16062647299186369205e-03;
          E[1][3]= 0.64818802653621866868e-02;
          E[1][4]=-0.12132386914873895089e-01;
          E[1][5]=-0.99709737725909584834e-02;
          E[1][6]=-0.70287093442894942752e-02;
          E[1][7]= 0.31243249755879001843e-01;
          E[1][8]= 0.31763603839792897936e-04;
          E[2][0]=-0.40875203230945019464e+00;
          E[2][1]= 0.28214948905763253599e+00;
          E[2][2]=-0.22612660499718519054e+00;
          E[2][3]= 0.13640993962985420478e+00;
          E[2][4]= 0.15888529591697266925e+00;
          E[2][5]=-0.11667863471317749710e+01;
          E[2][6]= 0.25224964119340060668e+00;
          E[2][7]= 0.10440940643938620983e+01;
          E[2][8]= 0.33914722176493324285e-03;
          E[3][0]=-0.29437531285359759661e+01;
          E[3][1]= 0.20017220470127690267e+01;
          E[3][2]=-0.15383035791443948798e+01;
          E[3][3]= 0.78114323215109899716e+00;
          E[3][4]= 0.13930345104184182146e+01;
          E[3][5]=-0.75958281612589849630e+01;
          E[3][6]= 0.18220129530415584951e+01;
          E[3][7]= 0.62663163493155487560e+01;
          E[3][8]= 0.54279630166374655267e-02;
          E[4][0]=-0.79572842006457093076e+01;
          E[4][1]= 0.53527892762707449170e+01;
          E[4][2]=-0.40049139768467199697e+01;
          E[4][3]= 0.18326058141135591515e+01;
          E[4][4]= 0.39753886181058367500e+01;
          E[4][5]=-0.19423696478604790213e+02;
          E[4][6]= 0.49362128400107292627e+01;
          E[4][7]= 0.15601708062381928560e+02;
          E[4][8]= 0.32142123424873719685e-01;
          E[5][0]=-0.78463118056075171475e+01;
          E[5][1]= 0.53580869574441241664e+01;
          E[5][2]=-0.41476905275607763365e+01;
          E[5][3]= 0.21275912797813913113e+01;
          E[5][4]= 0.37642416878253538582e+01;
          E[5][5]=-0.20329681631523484613e+02;
          E[5][6]= 0.48515418060343387549e+01;
          E[5][7]= 0.16604467346259915039e+02;
          E[5][8]= 0.84559690262225766975e-01;
          SM[0]= 0.00000000000000000000e+00;
          SM[1]= 0.10000000000000000000e+01;
          SM[2]= 0.17500000000000000000e+01;
          AM[0]= 0.58080578375796358720e+05;
          AM[1]=-0.33214989339522861968e+05;
          AM[2]= 0.28376088288312020853e+05;
          AM[3]=-0.27923430684614999462e+05;
          AM[4]= 0.29743005589491042677e+05;
          AM[5]=-0.15525927919158826444e+05;
          AM[6]=-0.27700591278076171875e+03;
          AM[7]= 0.73086943817138671875e+03;
          AM[8]= 0.00000000000000000000e+00;

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



} } 



