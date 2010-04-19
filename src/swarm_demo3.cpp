#include "swarm.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <memory>

#include<gsl/gsl_multifit.h>

const double DaysInYear = 365.2425;
const double TimeUnitToMetersPerSec = 29785.9163;
const double MetersPerSecToAuPerYear = 1./4740.57581;

double DrawUniform01() 
{ return static_cast<double>(rand())/static_cast<double>(RAND_MAX); }

double DrawStdNormal()
{
  double rn1 = DrawUniform01();
  double rn2 = DrawUniform01();
  return sqrt(-2.*log(rn1))*cos(2.*M_PI*(rn2));
}

double improve_mean_to_eccentric_annomaly_guess(const double e, const double M, const double x)
{
      double sx = sin(x);
      double cx = cos(x);
      double es = e*sx;
      double ec = e*cx;
      double f = x-es-M;
      double fp  = 1.-ec;
      double fpp = es;
      double fppp = ec;
      double dx = -f/fp;
      dx = -f/(fp+dx*fpp/2.);
      dx = -f/(fp+dx*fpp/2.+dx*dx*fppp/6.);
      return x+dx;
};

double mean_to_eccentric_annomaly(const double e,  double M)
{
  const int ORBEL_EHIE_NMAX = 3;

  int nper = static_cast<int>(M/(2.*M_PI));
  M = M - nper*2.*M_PI;
  if(M<0.)  M = M + 2.*M_PI;
  assert(M>=0.);
  if(M>M_PI) 
    { 
      M = 2.*M_PI - M;  
      double x = pow(6.*M,1./3.);
      for(int i=1;i<=ORBEL_EHIE_NMAX;++i)
	x = improve_mean_to_eccentric_annomaly_guess(e,M,x);
      x = 2.*M_PI-x;  
      return x;
    }
  else
    {
      double x = pow(6.*M,1./3.);
      for(int i=1;i<=ORBEL_EHIE_NMAX;++i)
	x = improve_mean_to_eccentric_annomaly_guess(e,M,x);
      return x;
    }
}

double DrawValue(const swarm::config& cfg, const std::string& name, const int bod, double min, double max)
{
  double val;
  std::stringstream s;

  s.str(""); s << name << '_' << bod << "_min";      
  if(cfg.count(s.str()))
    {
      double min_val = atof(cfg.at(s.str()).c_str());
      min = std::max(min,min_val);
    }
  s.str(""); s << name << '_' << bod << "_max";      
  if(cfg.count(s.str()))
    {
      double max_val = atof(cfg.at(s.str()).c_str());
      max = std::min(max,max_val);
    }

  s.str(""); s << name << '_' << bod;
  if(cfg.count(s.str()))
    {
      val = atof(cfg.at(s.str()).c_str());      
      if((val<min)||(val>max))
	{
	  std::cerr << "# " << s.str() << ": " << val << " ["<< min << ", " << max << "]\n";
	}
      assert(val>=min); 
      assert(val<=max);
      
      s.str(""); s << name << '_' << bod << "_sigma";      
      if(cfg.count(s.str()))
	{
	  double mean = val;
	  double sigma = atof(cfg.at(s.str()).c_str());
	  if(sigma)
	    {
	      do {
		val = mean + sigma*DrawStdNormal();
	      } while((val<min)||(val>max));
	    }
	  else
	    { val = mean; }
	}
    }
  else
    {
      val = min + (max-min)*DrawUniform01();
    }
  return val;
}

void calc_cartesian_for_ellipse(double& x,double& y, double & z, double &vx, double &vy, double &vz, const double a, const double e, const double i, const double O, const double w, const double M, const double GM)
{
  double cape = mean_to_eccentric_annomaly(e,M);
  double scap, ccap;
  sincos(cape,&scap,&ccap);
  double sqe = sqrt(1.-e*e);
  double sqgma = sqrt(GM*a);
  double xfac1 = a*(ccap-e);
  double xfac2 = a*sqe*scap;
  double ri = 1./(a*(1.-e*ccap));
  double vfac1 = -ri*sqgma * scap;
  double vfac2 =  ri*sqgma*sqe*ccap;

  double sw, cw, so, co, si, ci;
  sincos(w,&sw,&cw);
  sincos(O,&so,&co);
  sincos(i,&si,&ci);
  double d1[] = { cw*co-sw*so*ci, cw*so+sw*co*ci, sw*si};
  double d2[] = {-sw*co-cw*so*ci,-sw*so+cw*co*ci, cw*si};
  x  = d1[0]*xfac1+d2[0]*xfac2;
  y  = d1[1]*xfac1+d2[1]*xfac2;
  z  = d1[2]*xfac1+d2[2]*xfac2;
  vx = d1[0]*vfac1+d2[0]*vfac2;
  vy = d1[1]*vfac1+d2[1]*vfac2;
  vz = d1[2]*vfac1+d2[2]*vfac2;
}

void set_initial_conditions_for_demo(swarm::ensemble& ens, const swarm::config& cfg, const double t_first_obs) 
{

  bool use_jacobi = cfg.count("use_jacobi") ? atoi(cfg.at("use_jacobi").c_str()) : 0;

  std::cerr << "Set initial time for all systems = ";
  double time_init = cfg.count("time_init") ? atof(cfg.at("time_init").c_str()) : t_first_obs;
  std::cerr << time_init << ".\n";
  time_init *= 2.*M_PI/DaysInYear;
  ens.set_time_all(time_init);	

  for(unsigned int sys=0;sys<ens.nsys();++sys)
    {
      // set sun to unit mass and at origin
      float mass_star = cfg.count("mass_star") ? atof(cfg.at("mass_star").c_str()) : 1.;
      double x=0, y=0, z=0, vx=0, vy=0, vz=0;
      ens.set_body(sys, 0, mass_star, x, y, z, vx, vy, vz);

      double mass_enclosed = mass_star;
      for(unsigned int bod=1;bod<ens.nbod();++bod)
	{

	  double P = DrawValue(cfg,"P",bod,0.1,100000.);
	  double K = DrawValue(cfg,"K",bod,0.1,10000.);
	  double e = DrawValue(cfg,"ecc",bod,0.,1.);
	  double i = DrawValue(cfg,"inc",bod,-180.,180.);
	  i *= M_PI/180.;	  
          P /= DaysInYear;
          K *= MetersPerSecToAuPerYear;
	  double sini = sin(i);
          double mu = 0., a;
          for(unsigned int aa=0;aa<6;++aa)
            {
              a = pow((mass_enclosed+mu*mass_star)*P*P, 1.0/3.0);
              mu = K*P*sqrt(1.-e*e)/(2.*M_PI*a*sini);
            }
	  double mass_planet = mass_star*mu;
	  //	  float mass_planet = DrawValue(cfg,"mass",bod,0.,mass_star);
	  mass_enclosed += mass_planet;


	  double O = DrawValue(cfg,"node",bod,-720.,720.);
	  double w = DrawValue(cfg,"omega",bod,-720.,720.);
	  double M = DrawValue(cfg,"meananom",bod,-720.,720.);
	  //	  std::cout << "# Drawing sys= " << sys << " bod= " << bod << ' ' << mass_planet << "  " << a << ' ' << e << ' ' << i << ' ' << O << ' ' << w << ' ' << M << '\n';


	  O *= M_PI/180.;
	  w *= M_PI/180.;
	  M *= M_PI/180.;

	  double mass = use_jacobi ? mass_enclosed : mass_star+mass_planet;
	  calc_cartesian_for_ellipse(x,y,z,vx,vy,vz, a, e, i, O, w, M, mass);

	  double bx, by, bz, bvx, bvy, bvz;
	  ens.get_barycenter(sys,bx,by,bz,bvx,bvy,bvz,bod-1);
	  x  += bx;	  y  += by;	  z  += bz;
	  vx += bvx;	  vy += bvy;	  vz += bvz;
	  
	  // assign body a mass, position and velocity
	  ens.set_body(sys, bod, mass_planet, x, y, z, vx, vy, vz);
	}  // end loop over bodies

      // Shift into barycentric frame
      ens.get_barycenter(sys,x,y,z,vx,vy,vz);
      for(unsigned int bod=0;bod<ens.nbod();++bod)
	{
	  ens.set_body(sys, bod, ens.mass(sys,bod), 
		       ens.x(sys,bod)-x, ens.y(sys,bod)-y, ens.z(sys,bod)-z, 
		       ens.vx(sys,bod)-vx, ens.vy(sys,bod)-vy, ens.vz(sys,bod)-vz);	  
	}  // end loop over bodies
    } // end loop over systems
}

void calc_keplerian_for_cartesian( double& a,  double& e,  double& i,  double& O,  double& w,  double& M, const double x,const double y, const double z, const double vx, const double vy, const double vz, const double GM)
{
  const double TINY = 1.e-8;

  double h[] = {y*vz-z*vy, z*vx-x*vz, x*vy-y*vx};
  double h2 = h[0]*h[0]+h[1]*h[1]+h[2]*h[2];
  double hh = sqrt(h2);
  i = acos(h[2]/hh);
  double fac = sqrt(h[0]*h[0]+h[1]*h[1])/hh;
  double u;
  if(fac<TINY)
    {
      O = 0.;
      u = atan2(y,x);
      if(fabs(i-M_PI)<10.*TINY) u = -u;
    }
  else
    {
      O = atan2(h[0],-h[1]);
      u = atan2(z/sin(i),x*cos(O)+y*sin(O));
    }
  if(O<0.) O += 2.*M_PI;
  if(u<0.) u += 2.*M_PI;
  double r = sqrt(x*x+y*y+z*z);
  double energy = (vx*vx+vy*vy+vz*vz)*0.5-GM/r;
  
  if(fabs(energy*r/GM)<sqrt(TINY))
    { // Parabola
      a = 0.5*h2/GM;
      e = 1.;
      double ww = acos(2.*a/r-1.);
      if(vx*x+vy*y+vz*z<0.) w = 2.*M_PI-w;
      double tmpf = tan(0.5*w);
      M = tmpf*(1.+tmpf*tmpf/3.);
      w = u-ww;
      if(w<0.) w+= 2.*M_PI;
      w -= round(w/(2.*M_PI))*2.*M_PI;
    }
  else if (energy<0)
    { // Elipse
      a = -0.5*GM/energy;
      fac = 1.-h2/(GM*a);
      double ww, cape;
      if(fac>TINY)
	{
	  e = sqrt(fac);
	  double face = (a-r)/(a*e);
	  if(face>1.) cape = 0.;
	  else if (face>-1.) cape = acos(face);
	  else cape = M_PI;
	  
	  if(vx*x+vy*y+vz*z<0.) cape = 2.*M_PI-cape;
	  double cw = (cos(cape)-e)/(1.-e*cos(cape));
	  double sw = sqrt(1.-e*e)*sin(cape)/(1.-e*cos(cape));
	  ww = atan2(sw,cw);
	  if(ww<0.) ww += 2.*M_PI;
	}
      else
	{
	  e = 0.;
	  ww = u;
	  cape = u;
	}
      M = cape - e*sin(cape);
      w = u - ww;
      if(w<0.) w += 2.*M_PI;
      w -= round(w/(2.*M_PI))*2.*M_PI;
    }
  else if (energy>0)
    { // Hyperbola
      a = 0.5*GM/energy;
      fac = h2/(GM*a);
      double ww, capf;
      if(fac>TINY)
	{
	  e = sqrt(1.+fac);
	  double tmpf = (a+r)/(a*e);
	  capf = log(tmpf+sqrt(tmpf*tmpf-1.));
	  if(vx*x+vy*y+vz*z<0.) capf = -capf;
	  double cw = (e-cosh(capf))/(e*cosh(capf)-1.);
	  double sw = sqrt(e*e-1.)*sinh(capf)/(e*cosh(capf)-1.);
	  ww = atan2(sw,cw);
	  if(ww<0.) ww += 2.*M_PI;	  
	}
      else
	{
	  e = 1.;
	  double tmpf = 0.5*h2/GM;
	  ww = acos(2.*tmpf/r-1.);
	  if(vx*x+vy*y+vz*z<0.) ww = 2.*M_PI-ww;
	  tmpf = (a+r)/(a*e);
	  capf = log(tmpf+sqrt(tmpf*tmpf-1.));
	}
      M = e*sinh(capf)-capf;
      w = u - ww;
      if(w<0.) w+=2.*M_PI;
      w -= round(w/(2.*M_PI))*2.*M_PI;
    }
};

void print_selected_systems_for_demo(swarm::ensemble& ens)
{
  std::streamsize cout_precision_old = std::cout.precision();
  std::cout.precision(10);
  unsigned int nprint = std::max(10,ens.nsys());
  for(unsigned int systemid = 0; systemid< nprint; ++systemid)
    {
      std::cout << "sys= " << systemid << " time= " << ens.time(systemid) << "\n";
      double mass_enclosed = 0.;
      for(unsigned int bod=0;bod<ens.nbod();++bod)
	{
	  double mass = ens.mass(systemid,bod);
	  mass_enclosed += mass;

	  if(bod>0)
	    {
	      std::cout << "body= " << bod << ": ";
	      //	  std::cout << "pos= (" << ens.x(systemid, bod) << ", " <<  ens.y(systemid, bod) << ", " << ens.z(systemid, bod) << ") vel= (" << ens.vx(systemid, bod) << ", " <<  ens.vy(systemid, bod) << ", " << ens.vz(systemid, bod) << ").\n";

	      double bx, by, bz, bvx, bvy, bvz;
	      ens.get_barycenter(systemid,bx,by,bz,bvx,bvy,bvz,bod-1);
	      double x = ens.x(systemid,bod)-bx;
	      double y = ens.y(systemid,bod)-by;
	      double z = ens.z(systemid,bod)-bz;
	      double vx = ens.vx(systemid,bod)-bvx;
	      double vy = ens.vy(systemid,bod)-bvy;
	      double vz = ens.vz(systemid,bod)-bvz;

	      double a, e, i, O, w, M;
	      calc_keplerian_for_cartesian(a,e,i,O,w,M, x,y,z,vx,vy,vz, mass_enclosed);
	      O *= 180/M_PI;
	      w *= 180/M_PI;
	      M *= 180/M_PI;
	      std::cout << " a= " << a << " e= " << e << " i= " << i << " O= " << O << " w= " << w << " M= " << M << "\n";
	    }

	}
    }
  std::cout.precision(cout_precision_old);
}

class ObservationListClass
{
public:
  unsigned int mCount, mMaxObsId;
  std::vector<double> mTime, mRv, mErr;
  std::vector<int> mObsId;
  double mJitter;
  ObservationListClass() : mCount(0), mMaxObsId(0), mTime(0), mRv(0), mErr(0), mObsId(0) {};

  ObservationListClass(const char *Filename) : mCount(0), mMaxObsId(0), mTime(0), mRv(0), mErr(0), mObsId(0) 
  {
    Read(Filename);
  };

  void AddObservation(const double t, const double rv, const double err, const unsigned int obsid)
  {
    mTime.push_back(t);
    mRv.push_back(rv);
    mErr.push_back(err);
    assert(obsid<=mMaxObsId);
    mObsId.push_back(obsid);
    ++mCount;
  };

  void Read(const char *Filename)
  {
    //    std::cerr << "# Opening " << Filename << " to read observations.\n";
    std::ifstream ObsFile(Filename);
    if(!ObsFile.good()) { std::cerr << "Can't open " << Filename << " for reading observations\n"; exit(255); }
    for(mCount=0;ObsFile.good();)
      {
        char line[1024];
        ObsFile.getline(line,1023);
        if(!ObsFile.good()) { break; }
        std::istringstream LineStr;
        LineStr.str(line);
        double t, rv, err;
        unsigned int obs;
        LineStr >> t >> rv >> err >> obs;
        if(obs>mMaxObsId) mMaxObsId = obs;
        AddObservation(t,rv,err,obs);
        assert(obs>=0);
	std::cerr << "# Read observation: " << t << " " << rv << " " << err << " " << obs << "\n";
        
      } 
    ObsFile.close();
    std::cerr << "# Done reading " << mCount << " observations.\n";
  };

  double time(const int id) const  { return mTime[id];  };
  double data(const int id) const  { return mRv[id];    };
  double sigma(const int id) const { return mErr[id];   };
  int obsid(const int id) const    { return mObsId[id]; };
  int numobs() const               { return mCount; };
  int numobsids() const            { return mMaxObsId+1; };

};

class CalcChiSqForBestConstantsClass
{
  int num_obs, num_param;
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;
  gsl_multifit_linear_workspace * work;

public:
  CalcChiSqForBestConstantsClass(const int nobs, const int nparam) : num_obs(nobs), num_param(nparam)
  {
    allocate();
  };

  ~CalcChiSqForBestConstantsClass()
  {
    deallocate();
  };

  void allocate()
  {
    X = gsl_matrix_alloc (num_obs, num_param);
    y = gsl_vector_alloc (num_obs);
    w = gsl_vector_alloc (num_obs);
    c = gsl_vector_alloc (num_param);
    cov = gsl_matrix_alloc (num_param, num_param);
    work = gsl_multifit_linear_alloc (num_obs, num_param);
  };

  void deallocate()
  {
    gsl_matrix_free(X);
    gsl_vector_free(y);
    gsl_vector_free(w);
    gsl_vector_free(c);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free (work);
  };

  void set_predictors(const std::vector<int>& obsid) 
  {
    assert(obsid.size()==num_obs);
    for (int i = 0; i < num_obs; i++)
      {
	for(int j = 0; j < num_param; j++)
	  gsl_matrix_set(X, i,j, 0.);
	gsl_matrix_set(X,i,obsid[i],1.);
      }
  };

  void set_sigma(const std::vector<double>& sigma) 
  {
    assert(sigma.size()==num_obs);
    for (int i = 0; i < num_obs; i++)
      gsl_vector_set (w, i, 1.0/(sigma[i]*sigma[i]));
  };

  void set_data(const std::vector<double>& obs, const std::vector<double>& model)
  {
    assert(obs.size()==num_obs);
    assert(model.size()==num_obs);
    for (int i = 0; i < num_obs; i++)
      gsl_vector_set (y, i, obs[i]-model[i]);
  };

#if 0
  void setup(const ObservationListClass& data, const std::vector<double>& model)
  {
    set_predictors(data.obsid);
    set_sigma(data.sigma);
    set_data(data.data,model);
  };
#endif

  double operator()() 
  {
    double chisq;
    gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
    return chisq;
  };
  
  double param(const int i) const { return gsl_vector_get(c,i); };
  double covar(const int i, const int j) const { return gsl_matrix_get(cov,i,j); };

};

std::vector<double> CalcChiSq(swarm::cpu_ensemble& ens, swarm::gpu_ensemble& gpu_ens, std::auto_ptr<swarm::integrator> integ_gpu, const ObservationListClass& Obs)
{
  std::cerr << "Upload data to GPU.\n";
  gpu_ens.copy_from(ens);

  CalcChiSqForBestConstantsClass fit(Obs.numobs(),Obs.numobsids());
  fit.set_predictors(Obs.mObsId);
  fit.set_sigma(Obs.mErr);
  //	std::vector<std::vector<double> > StarVz(Obs.numobs(),std::vector<double>(nsystems));
  std::vector<std::vector<double> > StarVz(ens.nsys(),std::vector<double>(Obs.numobs()));
  for(int sysid=0;sysid<ens.nsys();++sysid)
    StarVz[sysid][0] = ens.vz(sysid,0);
  
  for(int obsid=1;obsid<Obs.numobs();++obsid)
    {
#if 0   // How should we deal with stop criteria here?
      std::cerr << "Set stop time for all systems.\n";
      double time_stop = Obs.time(obsid)*2.*M_PI/DaysInYear;
      ens.set_time_end_all(time_stop);
#endif
      
      double dT = (Obs.time(obsid)-Obs.time(obsid-1))*2.*M_PI/DaysInYear;
      
      // Perform the integration on gpu
      std::cerr << "Integrate ensemble on GPU.\n";
      integ_gpu->integrate(gpu_ens, dT);				
      std::cerr << "Download data to host.\n";
      ens.copy_from(gpu_ens);					
      std::cerr << "GPU integration complete.\n";
      // Save data for calculating chi^2
      for(int sysid=0;sysid<ens.nsys();++sysid)
	StarVz[sysid][obsid] = ens.vz(sysid,0);
    }
  
  // Check Energy conservation
  //  ens.calc_total_energy(&energy_final[0]);
  
  std::vector<double> ChiSq(ens.nsys());
  for(int sysid=0;sysid<ens.nsys();++sysid)
    {
      fit.set_data(Obs.mRv,StarVz[sysid]);
      ChiSq[sysid] = fit();
      //      double deltaE = (energy_final[sysid]-energy_init[sysid])/energy_init[sysid];
      std::cout << sysid << ": " << ChiSq[sysid] << "\n"; // " dE/E= " << deltaE << "\n";
    }
  
  return ChiSq;
}


void calc_semimajor_axes(const swarm::ensemble& ens, std::vector<std::vector<double> >& semimajor_axes)
{
  for(unsigned int systemid = 0; systemid< ens.nsys(); ++systemid)
    {
      //      std::cout << "sys= " << systemid << " time= " << ens.time(systemid) << "\n";
      double mass_enclosed = 0.;
      for(unsigned int bod=0;bod<ens.nbod();++bod)
	{
	  double mass = ens.mass(systemid,bod);
	  mass_enclosed += mass;

	  if(bod<=0) continue;
	  //	  std::cout << "body= " << bod << ": ";
	  //	  std::cout << "pos= (" << ens.x(systemid, bod) << ", " <<  ens.y(systemid, bod) << ", " << ens.z(systemid, bod) << ") vel= (" << ens.vx(systemid, bod) << ", " <<  ens.vy(systemid, bod) << ", " << ens.vz(systemid, bod) << ").\n";
	  
	  double bx, by, bz, bvx, bvy, bvz;
	  ens.get_barycenter(systemid,bx,by,bz,bvx,bvy,bvz,bod-1);
	  double x = ens.x(systemid,bod)-bx;
	  double y = ens.y(systemid,bod)-by;
	  double z = ens.z(systemid,bod)-bz;
	  double vx = ens.vx(systemid,bod)-bvx;
	  double vy = ens.vy(systemid,bod)-bvy;
	  double vz = ens.vz(systemid,bod)-bvz;
	  
	  double a, e, i, O, w, M;
	  calc_keplerian_for_cartesian(a,e,i,O,w,M, x,y,z,vx,vy,vz, mass_enclosed);
	  O *= 180/M_PI;
	  w *= 180/M_PI;
	  M *= 180/M_PI;
	  semimajor_axes[systemid][bod] = a;
	  //	  std::cout << " a= " << a << " e= " << e << " i= " << i << " O= " << O << " w= " << w << " M= " << M << "\n";
	}
    }
}

std::vector<int> choose_systems_to_halt(const swarm::ensemble& ens, const std::vector<std::vector<double> >& a_init)
{
  std::vector<std::vector<double> > a_final(ens.nsys(),std::vector<double>(ens.nbod()-1));
  calc_semimajor_axes(ens,a_final);
  std::vector<int> halt(ens.nsys(),0);
  for(int sysid=0;sysid<ens.nsys();++sysid)
    {
      int bad = 0;
      for(int plid=0;plid<ens.nbod();++plid)
	{
	  if(fabs(a_final[sysid][plid]-a_init[sysid][plid])>0.5*a_init[sysid][plid])
	    ++bad;
	}
      if(bad!=0) halt[sysid] = 1;
    }
  return halt;
}

int main(int argc, const char **argv)
{
  using namespace swarm;

  if(argc != 2)
    {
      std::cerr << "Usage: " << argv[0] << " <integrator.cfg>\n";
      return -1;
    }
  std::string icfgfn = argv[1];

  config cfg;
  load_config(cfg,icfgfn);

  if(!cfg.count("observation_file")) 
    {
      std::cerr << "# Must specify observation_file in configuration file.\n";
      return -1;
    }
  ObservationListClass Obs(cfg.at("observation_file").c_str());
  assert(Obs.numobs()>2);

  std::cerr << "# Initialize ensemble on host to be used with GPU integration.\n";
  unsigned int nsystems = cfg.count("num_systems") ? atoi(cfg.at("num_systems").c_str()) : 1024;
  unsigned int nbodyspersystem = cfg.count("num_bodies") ? atoi(cfg.at("num_bodies").c_str()) : 3;
  cpu_ensemble ens(nsystems, nbodyspersystem);
  gpu_ensemble gpu_ens(nsystems, nbodyspersystem);

	std:: cerr << "Initialize the GPU integrator\n";
	std::auto_ptr<integrator> integ_gpu(integrator::create(cfg));

	std::cerr << "Set initial conditions on CPU.\n";
	set_initial_conditions_for_demo(ens,cfg,Obs.time(0));

	// Print initial conditions on CPU for use w/ GPU
	//	std::cerr << "Print selected initial conditions for upload to GPU.\n";
	//	print_selected_systems_for_demo(ens);
	
	std::vector<double> energy_init(nsystems), energy_final(nsystems);
	ens.calc_total_energy(&energy_init[0]);

	std::vector<std::vector<double> > a_init(ens.nsys(),std::vector<double>(ens.nbod()-1));
	calc_semimajor_axes(ens,a_init);
	std::vector<double> ChiSq(ens.nsys()); //= CalcChiSq(ens,gpu_ens,integ_gpu,Obs);
  gpu_ens.copy_from(ens);

	// Now do stability test
	int num_integrations = 3;
	for(int intid=0;intid<num_integrations;++intid)
	  {
	    double time_stop = 100.*2.*M_PI*pow(2,intid); // +Obs.time(0)/DaysInYear;
	    double dT = time_stop-ens.time(0);
	    std::cerr << "Integrate ensemble on GPU from t= " << ens.time(0) << " to " << time_stop << ".\n";
	    integ_gpu->integrate(gpu_ens, dT);				
	    std::cerr << "Download data from GPU\n";
	    ens.copy_from(gpu_ens);					
	    std::cerr << "Print selected systems\n";
	    print_selected_systems_for_demo(ens);
	    std::cerr << "Choose Systems to halt\n";
	    std::vector<int> halt_flag = choose_systems_to_halt(ens,a_init);
	    std::cerr << "Set Inactive Flag\n";
	    ens.set_inactive(halt_flag);
	    std::cerr << "Pack CPU ensemble\n";
	    int nactive = ens.pack();
	  }


	// Check Energy conservation
	ens.calc_total_energy(&energy_final[0]);
	double max_deltaE = 0;
	for(int sysid=0;sysid<ens.nsys();++sysid)
	  {
	    double deltaE = (energy_final[ens.index_of_system(sysid)]-energy_init[sysid])/energy_init[sysid];
	    if(fabs(deltaE)>max_deltaE)
	      { max_deltaE = fabs(deltaE); }
	    std::cout << sysid << ": " << ChiSq[ens.index_of_system(sysid)] << " dE/E= " << deltaE << "\n";
	  }
	std::cerr << "# Max dE/E= " << max_deltaE << "\n";

	// Print results
	std::cerr << "Print selected results from GPU's calculation.\n";
	print_selected_systems_for_demo(ens);

	// both the integrator & the ensembles are automatically deallocated on exit
	// so there's nothing special we have to do here.
	return 0;
}

