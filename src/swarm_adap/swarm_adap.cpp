/*
    "swarm_adap" is a program that uses the Swarm-NG tools for modeling an ensemble of
    small N systems using the hermite_adap_gpu integrator.
    Copyright (C) 2010  Swarm-NG Development Group

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact aaron.boley at the domain gmail.com if you have questions regarding this
    software.

    SETUP: Initial conditions are generated using scatter.py.  This should be
    provided with the distribution, and resides in the "ic" directory. In the
    configuration file, which the code expects to be called integration-scatter.cfg, 
    make sure that h, pre, and stepfac are defined.  Be careful with stepfac, it 
    should be set to about 0.03 based on some minimal testing.  Feel free to experiment. 

*/
#include "swarm.h"
#include <memory>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>
#include <iomanip>
using namespace std;
using namespace swarm;

#include "swarm_adap/swarm_adap.h"

int main()
{
        cout.setf(ios::scientific,ios::floatfield);
	// load the ensemble
	cpu_ensemble ens;
	load_ensemble("../scripts/ic/data", ens);
        unsigned int nSystems=ens.nsys(); 
	vector<double> E(nSystems,0.);

        cout<<" nSystems "<<nSystems<<"\n";

	// set up the integrator and integrator config (TODO: load from config file)
	config cfg;
	load_config(cfg, "integrator.cfg");
        check_cfg_input(ens,cfg);
        init(cfg);
	std::auto_ptr<integrator> integ(integrator::create(cfg));

        //get Observing Times
        vector<real>   ObsTimes=getObsTimes("../scripts/ic/");
        unsigned int nObs=ObsTimes.size();

        //open and clear files
        for (unsigned int i=0;i<nSystems;++i)
         {
          ofstream thisFile;
          stringstream buff;

          buff.str("");
          buff<<OUTPUT_DIRECTORY<<"/systemInfo."<<setfill('0')<<setw(5)<<i;
          thisFile.open(buff.str().c_str(),ios::trunc);
          assert(thisFile.good());
          thisFile.close();
         }


        //
        // Set timing and disable logging -- will do it by hand below
        //

        unsigned int observation=0;
        ens.set_time_all(ObsTimes[observation]); // initial time
        ens.set_time_end_all(ObsTimes[nObs-1]);  // max integration time
        ens.set_time_output_all(1, 1.01*ObsTimes[nObs-1]);

        gpu_ensemble gpu_ens(ens);				// upload to GPU

        //
        //Now integrate, but stop at each observation time to check progress and log data.
        //

 
        real startTime=ObsTimes[observation];
        while(observation++<nObs-1)
         {
           real dT=ObsTimes[observation]-startTime;
        // check progress by writing to stdout
           cout<<" Working on observation "<<observation<<" of "<<nObs-1<<' ' <<"with time interval "<<dT/yrToCodeTime<<" in yr \n";

           integ->integrate(gpu_ens, dT);				// integrate
           cudaThreadSynchronize();
           ens.copy_from(gpu_ens);					// download to host
           startTime=ObsTimes[observation];

	// store output
           for (unsigned int i=0;i<nSystems;++i) 
            {
              ofstream thisFile;
              stringstream buff;
  
              buff.str("");
              buff<<OUTPUT_DIRECTORY<<"/systemInfo."<<setfill('0')<<setw(5)<<i;
              thisFile.open(buff.str().c_str(),ios::app);
              assert(thisFile.good());
              thisFile<<scientific<<setprecision(9)<<ens.time(i)/yrToCodeTime<<' ';
              double Enew = calc_system_energy(ens,i);
              for(int bod = 0; bod != ens.nbod(); bod++)
              {
                float m; double x[3], v[3];
                ens.get_body(i, bod, m, x[0], x[1], x[2], v[0], v[1], v[2]);
                thisFile<<scientific<<setprecision(9)<<m<<' '<<x[0]<<' '<<x[1]<<' '<<x[2]<<' '<<v[0]/kmpsToCodeVel<<' '<<v[1]/kmpsToCodeVel<<' '<<v[2]/kmpsToCodeVel<<' ';

              }
              thisFile<<Enew<<"\n";
              thisFile.close();

            }
 

          } // end while loop

	return 0;
}
