/*
    Header file for the program swarm_adap.
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

    Contact aaron.boley at the domain gmail.com if you have questions regarding
    this software.
*/

#ifndef __SWARM_SCATTER_DEMO_H__
#define __SWARM_SCATTER_DEMO_H__

#ifndef __SWARM_USER_H__
#define __SWARM_USER_H__
#include "user.h"
#endif

///////////////////////////////////////////
// The following can be overridden in user.h
// Or on command line.  These are defaults.
#ifndef OUTPUT_PREFIX 
#define OUTPUT_PREFIX "./scatter_output/systemInfo"
#endif

#ifndef OBSERVATION_FILE
#define OBSERVATION_FILE "../scripts/ic/observeTimes.dat"
#endif

//include directory with prefix
#ifndef IC_FILE_PREFIX
#define IC_FILE_PREFIX "../scripts/ic/data"
#endif

#ifndef CFG_FILE
#define CFG_FILE "integrator.cfg"
#endif
///////////////////////////////////////////

#define real double

#define GCGS     6.67e-8
#define AUCGS    1.496e13
#define yrToCodeTime (2.*M_PI)
#define secondsToCodeTime (2.*M_PI)/(365.25*24.*3600.)
#define kmpsToCodeVel (1./(AUCGS*1e-5*secondsToCodeTime))

// command line options
bool readCommandLine(int argc,char *argv[], string& cfgFileName, string& obsFileName, string& icsFileName, string& outFileName){
     if(argc>1){
       if((argc-1)%2==0){ // make sure there is a key-value pair.
         for(int i=1; i<argc; i+=2) {
           stringstream key,value;
           key<<argv[i];
           value<<argv[i+1];
           key.ignore(1024,' ');
           if(key.str()=="-cfg") {
		value>>cfgFileName;
           } else if(key.str()=="-obs") {
		value>>obsFileName;
           } else if(key.str()=="-ics") {
		value>>icsFileName;
           } else if(key.str()=="-out") {
		value>>outFileName;
           } else {
		cout<<" Usage: ./swarm_scatter_demo [ -cfg <configuration file> ], [ -obs <observation times file prefix> ], [ -ics <ic file prefix> ], [ -out <out file prefix> ]\n";
		cout<<" Example: ./swarm_scatter_demo -cfg integrator.cfg -obs ../scripts/ic/observeTimes.dat -ics ../scripts/ic/data -out scatter_output/systemInfo\n\n" ;
             return false;
           }

         }
       }
       else {
                cout<<" Usage: ./swarm_scatter_demo [ -cfg <configuration file> ], [ -obs <observation times file prefix> ], [ -ics <ic file prefix> ], [ -out <out file prefix> ] \n";
                cout<<" Example: ./swarm_scatter_demo -cfg integrator.cfg -obs ../scripts/ic/observeTimes.dat -ics ../scripts/ic/data -out scatter_output/systemInfo\n\n"              ;
		return false;
       }
     }
     else {
       cout << " Using default values (set at compile time). "<<endl;
     }

     cout<<" Using the following configuration file:" <<cfgFileName<<"\n";

     return true;
}
// let's check the configuration file
void check_cfg_input(config &cfg, string& obsFileName, string& icsFilePrefix, string& outFilePrefix)
{
	obsFileName = cfg.count("obsFileName") ? cfg["obsFileName"] : obsFileName;
	icsFilePrefix = cfg.count("icsFilePrefix") ? cfg["icsFilePrefix"] : icsFilePrefix;
	outFilePrefix = cfg.count("outFilePrefix") ? cfg["outFilePrefix"] : outFilePrefix; 
	std::string runon = cfg.count("runon") ? cfg["runon"] : "gpu";
        bool ongpu;
        if     (runon == "gpu") { ongpu = true; }
        else if(runon == "cpu") { ongpu = false; }
        else { ERROR("The 'runon' configuration file parameter must be one of 'gpu' or 'cpu'"); }
        std::string integrator_name = cfg["integrator"];
        std::cout << " Integrator: " << integrator_name << ", executing on the " << (ongpu ? "GPU" : "CPU") << "\n";
        if(integrator_name == "gpu_hermite_adap")
        { std::cout<<" Using minimum time step h = "<<cfg["h"]<<" and stepfac = "<<cfg["stepfac"]<<'\n'; } 
        else { std::cout<<" Using time step h = "<<cfg["h"]<<'\n'; }
        std::cout<<" Using rmax = "<<cfg["rmax"]<<'\n';
	std::cout<<" Using obsFileName = "<<obsFileName<<"\n";
	std::cout<<" Getting ICs from "<<icsFilePrefix<<"\n";
	std::cout<<" Writing output to "<<outFilePrefix<<"\n";

}

// compute the energy of a given system
double calc_system_energy(const cpu_ensemble &ens, const int sys)
{
  double E = 0.;
      for(int bod1 = 0; bod1 != ens.nbod(); bod1++)
	{
	  float m1; double x1[3], v1[3];
	  ens.get_body(sys, bod1, m1, x1[0], x1[1], x1[2], v1[0], v1[1], v1[2]);
	  E += 0.5*m1*(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
	  
	  for(int bod2 = 0; bod2 < bod1; bod2++)
	    {
	      float m2; double x2[3], v2[3];
	      ens.get_body(sys, bod2, m2, x2[0], x2[1], x2[2], v2[0], v2[1], v2[2]);
	      double dist = sqrt((x2[0]-x1[0])*(x2[0]-x1[0])+(x2[1]-x1[1])*(x2[1]-x1[1])+(x2[2]-x1[2])*(x2[2]-x1[2]));
	      
	      E -= m1*m2/dist;
	    }
	}
	return E;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Read observing file
vector<real> getObsTimes(string fileName)
 {
     ifstream ThisFile;
     stringstream buffer;
     string line;
     vector<real> ObsTimes;

     buffer.str("");//clear buffer
     ThisFile.open(fileName.c_str());
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
     ThisFile.close();
     return ObsTimes;

 }
////////////////////////////////////////////////////////////////////////
#endif
