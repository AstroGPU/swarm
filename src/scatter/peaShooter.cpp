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

#define OUTPUT_DIRECTORY "output_files"
#define real double

#define GCGS     6.67e-8
#define AUCGS    1.496e13
#define yrToCodeTime (2.*M_PI)
#define secondsToCodeTime (2.*M_PI)/(365.25*24.*3600.)
#define kmpsToCodeVel (1./(AUCGS*1e-5*secondsToCodeTime))

// compute the energies of each system
// NOTE: This is different from the one used for Euler integrator
double calc_system_energy(const cpu_ensemble &ens, const int sys)
{
  double E = 0.;
  
/*  for(int sys = 0; sys != ens.nsys(); sys++)
    {*/
//      E = 0.;
      for(int bod1 = 0; bod1 != ens.nbod(); bod1++)
	{
	  float m1; double x1[3], v1[3];
	  ens.get_body(sys, bod1, m1, x1[0], x1[1], x1[2], v1[0], v1[1], v1[2]);
	  E += 0.5*m1*(v1[0]*v1[0]+v1[1]*v1[1]+v1[3]*v1[3]);
	  
	  for(int bod2 = 0; bod2 < bod1; bod2++)
	    {
	      float m2; double x2[3], v2[3];
	      ens.get_body(sys, bod2, m2, x2[0], x2[1], x2[2], v2[0], v2[1], v2[2]);
	      double dist = sqrt((x2[0]-x1[0])*(x2[0]-x1[0])+(x2[1]-x1[1])*(x2[1]-x1[1])+(x2[2]-x1[2])*(x2[2]-x1[2]));
	      
	      E -= m1*m2/dist;
	    }
	}
	return E;
/*    }*/
}

// just a simple dump to stdout of one system
void write_output(const cpu_ensemble &ens, const int sys, double& Eold)
{
  double Enew = calc_system_energy(ens,sys);
  for(int bod = 0; bod != ens.nbod(); bod++)
    {
      float m; double x[3], v[3];
      ens.get_body(sys, bod, m, x[0], x[1], x[2], v[0], v[1], v[2]);
      
      printf("%5d %5d  T=%f  m=%f  pos=(% 9.5f % 9.5f % 9.5f)  vel=(% 9.5f % 9.5f % 9.5f)  E=%g", sys, bod, ens.time(sys), m, x[0], x[1], x[2], v[0], v[1], v[2],Enew);
      if(Eold!=0)
	printf("  dE/E=%g\n",(Enew-Eold)/Eold);
      else
	printf("\n");
    }
  Eold = Enew;
}


////////////////////////////////////////////////////////////////////////
// write log function
int writeLogToFile(string prefix, string log, unsigned int i)
{
        ofstream thisFile;
        stringstream buff;

        buff<<OUTPUT_DIRECTORY<<"/"<<prefix<<setfill('0')<<setw(4)<<i;
        thisFile.open(buff.str().c_str());
        assert(thisFile.good());
        thisFile<<log<<"\n";
        thisFile.close();

        return 1;
}
////////////////////////////////////////////////////////////////////////



//ACB///////////////////////////////////////////////////////////////////
// Read observing file
vector<real> getObsTimes()
 {
     ifstream ThisFile;
     stringstream buffer;
     string line;
     vector<real> ObsTimes;

     buffer.str("");//clear buffer
     //buffer<<IC_DIRECTORY<<"/"<<OBSERVE_TIME_FILE;
     buffer<<"../ic"<<"/"<<"observeTimes.dat";
     ThisFile.open(buffer.str().c_str());
     assert(ThisFile.good());

     while(!ThisFile.eof()) 
      {
        stringstream value;
        float floatTemp;
        getline(ThisFile,line);
        value<<line;
        if (ThisFile.eof())continue;
        value>>floatTemp;
        ObsTimes.push_back(floatTemp*yrToCodeTime);
      }

#if VERBOSE_OUTPUT>0
     for(unsigned int i=0;i<ObsTimes.size();++i)
      {
        cout<<ObsTimes[i]<<" Observation Time (code)\n";
      }
#endif
     //if(ObsTimes.size()<nObs){
     //   cout<<" Incompatible observation sizes.  Expected "<<nObs<<", but got "<<ObsTimes.size()<<endl;
     //   cout<<" Expected number should be <= actual. Recovering by increasing nObs."<<endl;
      //  nObs=ObsTimes.size();
     //}
     ThisFile.close();
     return ObsTimes;

 }
////////////////////////////////////////////////////////////////////////


int main()
{
        cout.setf(ios::scientific,ios::floatfield);
	// load the ensemble
	cpu_ensemble ens;
	load_ensemble("../ic/data", ens);
        unsigned int nSystems=ens.nsys(); 
	vector<double> E(nSystems,0.);

        cout<<" nSystems "<<nSystems<<"\n";
	//printf("Initial conditions...\n");
	//write_output(ens, 0, E[0]);
	//write_output(ens, 1, E[1]);
	//write_output(ens, 2, E[2]);

	// set up the integrator and integrator config (TODO: load from config file)
	config cfg;
	load_config(cfg, "integrator-acboley.cfg");
	std::auto_ptr<integrator> integ(integrator::create(cfg));
	std::string runon = cfg.count("runon") ? cfg["runon"] : "gpu";
	bool ongpu;
	     if(runon == "gpu") { ongpu = true; }
	else if(runon == "cpu") { ongpu = false; }
	else { ERROR("The 'runon' configuration file parameter must be one of 'gpu' or 'cpu'"); }
	std::cerr << "Integrator: " << cfg["integrator"] << ", executing on the " << (ongpu ? "GPU" : "CPU") << "\n";

        //ACB
        //get Observing Times
        vector<real>   ObsTimes=getObsTimes();
        unsigned int nObs=ObsTimes.size();
        //ACB

	/* Replace GPU call here with CPU call below comment
	   gpu_ensemble gpu_ens(ens);
	   integ->integrate(gpu_ens, dT);
	   ens.copy_from(gpu_ens);
	*/

        //open and clear files
        for (unsigned int i=0;i<nSystems;++i)
         {
          ofstream thisFile;
          stringstream buff;

          buff.str("");
          buff<<OUTPUT_DIRECTORY<<"/systemInfo."<<setfill('0')<<setw(4)<<i;
          thisFile.open(buff.str().c_str(),ios::trunc);
          assert(thisFile.good());
          thisFile.close();
         }


        //
        //ACB
        //Now integrate, but stop at each observation time to check progress and log data.
        //

        gpu_ensemble gpu_ens(ens);				// upload to GPU
        unsigned int observation=0;
        real startTime=ObsTimes[observation];
        while(observation++<nObs-1)
         {
           real dT=ObsTimes[observation]-startTime;
           cout<<" Working on observation and dt "<<observation<<" of "<<nObs-1<<' ' <<dT<<"\n";
           integ->integrate(gpu_ens, dT);				// integrate
           cudaThreadSynchronize();
           ens.copy_from(gpu_ens);					// download to host
           //cudaThreadSynchronize(); 
           startTime=ObsTimes[observation];

	// store output
           for (unsigned int i=0;i<nSystems;++i) 
            {
              ofstream thisFile;
              stringstream buff;
  
              buff.str("");
              buff<<OUTPUT_DIRECTORY<<"/systemInfo."<<setfill('0')<<setw(4)<<i;
              thisFile.open(buff.str().c_str(),ios::app);
              assert(thisFile.good());
              thisFile<<scientific<<setprecision(9)<<ens.time(i)<<' ';
              double Enew = calc_system_energy(ens,i);
              //logTheseData(ens, const int sys, string &log);
              for(int bod = 0; bod != ens.nbod(); bod++)
              {
                float m; double x[3], v[3];
                ens.get_body(i, bod, m, x[0], x[1], x[2], v[0], v[1], v[2]);
                thisFile<<scientific<<setprecision(9)<<m<<' '<<x[0]<<' '<<x[1]<<' '<<x[2]<<' '<<v[0]<<' '<<v[1]<<' '<<v[2]<<' ';

              }
              thisFile<<Enew<<"\n";
              thisFile.close();

            }
 

          } // end while loop

////////////////////////////////////////////////////////////////////////
// Output section. Write log streams to file.
//    for (unsigned int i=0;i<nSystems;++i) {
//        assert(writeLogToFile("logFile.",logData[i].str(),i));
//    }

	// both the integrator & the ensembles are automatically deallocated on exit
	// so there's nothing special we have to do here.
	return 0;
}
