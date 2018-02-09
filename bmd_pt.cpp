#include "bmd_pt.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
using namespace std;

const double ln2pi=1.83787706641;
const double ln2=0.69314718056;


inline double normsquared(const double x,const double y){
  return x*x+y*y;
}

double percentile(vector<double> x,double p){
	sort(x.begin(),x.end());
	return x[p*x.size()];
}

inline int integraltable::get_a_ind(const double a){
 return 1000.0/12.0*(log(a)+6.0)-1;
}

inline double integraltable::get_a(const int a_ind){
 return exp(((a_ind+1)*12.0/1000.0)-6.0);
}

inline int integraltable::get_d_ind(const double d){
   return (d-0.1)/d_step_size-1;
}

inline double integraltable::get_d(const int d_ind){
 return (d_ind+1)*d_step_size+0.1;
}

inline double integraltable::operator()(const double a,const double d){
  int d_ind=get_d_ind(d), a_ind;
  double I1,I2,I3,I4;
  double a_diffup,a_diffdown,d_diffup,d_diffdown;
  if(a<a_lower[d_ind])
    return 0.25/a-d*log(a)-d*ln2;
  if (a>a_upper[d_ind])
    return -ln2-0.5*(d+1.0)*log(a)+log(tgamma(0.5*(d+1.0))+0.25*tgamma(0.5*(d+3.0))/a);
  else{
    a_ind=get_a_ind(a);
    I1=table[a_ind*Nd+d_ind];
    I2=table[a_ind*Nd+d_ind+1];
    I3=table[(a_ind+1)*Nd+d_ind];
    I4=table[(a_ind+1)*Nd+d_ind+1];
    a_diffup=get_a(a_ind+1)-a;
    d_diffup=get_d(d_ind+1)-d;
    a_diffdown=a-get_a(a_ind);
    d_diffdown=d-get_d(d_ind);
    return (d_diffdown*(a_diffdown*I4+a_diffup*I3) + d_diffup*(a_diffdown*I2+a_diffup*I1))
            /(a_diffup+a_diffdown)/(d_diffup+d_diffdown);
  }
}

integraltable::integraltable(char* filename){
  ifstream input_table(filename,ios::in);
  stringstream ss;
  string buff;
  d_step_size=0.0049;
  if(input_table.is_open()){
    getline(input_table, buff);
    ss << buff;
    ss >> Na;
    ss >> Nd;
    table=new double[Na*Nd];
    a_lower=new double[Nd];
    a_upper=new double[Nd];
    for(int row = 0; row < Na;row++) {
      getline(input_table,buff);
      ss.clear();
      ss << buff;
      for(int col = 0; col < Nd;col++)
        ss>>table[row*Nd+col];
    }
    getline(input_table,buff);
    ss.clear();
    ss << buff;
    for (int i=0; i<Nd; i++)
      ss>>a_lower[i];
    getline(input_table,buff);
    ss.clear();
    ss << buff;
    for (int i=0; i<Nd; i++)
       ss>>a_upper[i];
    input_table.close();
  }
}

integraltable::~integraltable(){
  delete[] table;
  delete[] a_lower;
  delete[] a_upper;
}

params::params(){
  //initialize velocity distribution parameters
  //d1=4.4;
  //sigma0=0.0003;
  //sigma1=0.03;
  // start from random parameters
  d1=uniform_real_distribution<double>{1.1,5.0}(engine);
  sigma0=uniform_real_distribution<double>{0.0001,0.005}(engine);
  sigma1=uniform_real_distribution<double>{0.005,0.1}(engine);
  sigmaz=0.015;
  sigmax=0.02;
  const_term_down=(0.5-0.5*d0)*ln2-(d0+1)*log(sigma0)-lgamma(0.5*(d0+1))+2.0*(d0+1.0)*log(sigmaz);
  const_term_up=(0.5-0.5*d1)*ln2-(d1+1)*log(sigma1)-lgamma(0.5*(d1+1))+2.0*(d1+1.0)*log(sigmaz);
  prefactor=ln2pi+2.0*log(sigmaz);
}

params::params(const double d1_val, const double sigma0_val, const double sigma1_val,
               const double sigmax_val, const double sigmaz_val){
  d1=d1_val;
  sigma0=sigma0_val;
  sigma1=sigma1_val;
  sigmax=sigmax_val;
  sigmaz=sigmaz_val;
  const_term_down=(1.0-0.5*(d0+1))*ln2-(d0+1)*log(sigma0)-lgamma(0.5*(d0+1))+2.0*(d0+1.0)*log(sigmaz);
  const_term_up=(1.0-0.5*(d1+1))*ln2-(d1+1)*log(sigma1)-lgamma(0.5*(d1+1))+2.0*(d1+1.0)*log(sigmaz);
  prefactor=ln2pi+2.0*log(sigmaz);
}

params::params(char* filename){
  ifstream input(filename,ios::in);
  input>>d1>>sigma0>>sigma1>>sigmaz>>sigmax;
  const_term_down=(1.0-0.5*(d0+1))*ln2-(d0+1)*log(sigma0)-lgamma(0.5*(d0+1))+2.0*(d0+1.0)*log(sigmaz);
  const_term_up=(1.0-0.5*(d1+1))*ln2-(d1+1)*log(sigma1)-lgamma(0.5*(d1+1))+2.0*(d1+1.0)*log(sigmaz);
  prefactor=ln2pi+2.0*log(sigmaz);
  input.close();
}

void params::set_sigmax(const double sigmax_val){
  sigmax=sigmax_val;
}

void params::set_sigmaz(const double sigmaz_val){
  sigmaz=sigmaz_val;
  const_term_down=(1.0-0.5*(d0+1))*ln2-(d0+1)*log(sigma0)-lgamma(0.5*(d0+1))+2.0*(d0+1.0)*log(sigmaz);
  const_term_up=(1.0-0.5*(d1+1))*ln2-(d1+1)*log(sigma1)-lgamma(0.5*(d1+1))+2.0*(d1+1.0)*log(sigmaz);
  prefactor=ln2pi+2.0*log(sigmaz);
}

void params::set_d_sigma_up(const double d_val, const double sigma_val){
  d1=d_val;
  sigma1=sigma_val;
  const_term_up=(1.0-0.5*(d1+1))*ln2-(d1+1)*log(sigma1)-lgamma(0.5*(d1+1))+2.0*(d1+1.0)*log(sigmaz);
}

void params::set_sigma_down(const double sigma_val){
  sigma0=sigma_val;
  const_term_down=(1.0-0.5*(d0+1))*ln2-(d0+1)*log(sigma0)-lgamma(0.5*(d0+1))+2.0*(d0+1.0)*log(sigmaz);
}

void params::write(){
  cout<<sigmaz<<"\t"<<sigmax<<"\t"<<d0<<"\t"<<sigma0<<"\t"<<d1<<"\t"<<sigma1<<endl;
}

void params::write(char* filename){
  ofstream output(filename,ios::out | ios::app);
  output<<sigmaz<<"\t"<<sigmax<<"\t"<<d0<<"\t"<<sigma0<<"\t"<<d1<<"\t"<<sigma1<<endl;
  output.close();
}

void data::kalman_filter(const params& theta){
  double sigzsq=theta.sigmaz*theta.sigmaz;
  double sigxsq=theta.sigmax*theta.sigmax;
  double P=0.5*(sqrt(sigzsq*sigzsq+4.0*sigzsq*sigxsq)-sigzsq);
  double K=(P+sigzsq)/(P+sigzsq+sigxsq);
  for(int dim=0;dim<2;dim++){
    z[dim][0]=x[dim][0];
    for(int i=1;i<T;i++)
      z[dim][i]=z[dim][i-1]+K*(x[dim][i]-z[dim][i-1]);
  }
  sum_dz_squared=0.0;
  for(int i=1;i<T;i++)
    sum_dz_squared+=normsquared(z[0][i]-z[0][i-1],z[1][i]-z[1][i-1]);
}

changepoints::changepoints(int T){
  t01.push_back(0);
  t01.push_back(T-1);
  t10.push_back(1);
  t10.push_back(T);
  N=2; n=2;
}

vector<bool> changepoints::get_vec(){
  vector<bool> C;
  C.insert(C.begin(),t10[0],1);
  for (int j=1;j<n;j++){
    C.insert(C.end(),t01[j]-t10[j-1],0);
    C.insert(C.end(),t10[j]-t01[j],1);
  }
  return C;
}

data::data(char* filename){
  ifstream input(filename,ios::in);
  input>>T;
  x[0]=new double[T];
  x[1]=new double[T];
  if(input.is_open())
  	for (int i=0;i<T;i++)
      input>>x[0][i]>>x[1][i];
  z[0]=new double[T];
  z[1]=new double[T];
  input.close();
}

data::~data(){
  delete x[0];
  delete x[1];
  delete z[0];
  delete z[1];
}

MCMC_chain::MCMC_chain(double beta_val,data* dat_ptr, params& theta_val, changepoints& C_val, integraltable* table_ptr):
beta(beta_val), dat(dat_ptr), theta(theta_val), C(C_val), table(table_ptr), get_step_type{1.0,1.0,1.0,1.0,2.0,2.0}{
  logprior=get_logprior();
  loglik=get_loglik();
  logposterior=logprior+loglik;
  logpostbest=logposterior;
  nsteps=0;
  swaps_acc=0;
  swaps_tried=0;
}

changepoints::changepoints(char* filename, int T){
  ifstream input(filename,ios::in);
  int C,Cnext=0;
  n=0,N=0;
  for(int t=0;t<T;t++){
    C=Cnext,input>>Cnext;
    if(C==0 && Cnext==1)
      t01.push_back(t),n++;
    if(C==1 && Cnext==0)
      t10.push_back(t);
    N+=Cnext;
  }
  t10.push_back(T);
  input.close();
}

changepoints::changepoints(data& dat){
  vector<double> dx_squared(dat.T,0.0);
  double threshold;
  n=1,N=1;
  dx_squared[0]=normsquared(dat.x[0][0],dat.x[1][0]);
  for(int t=1;t<dat.T;t++)
    dx_squared[t]=normsquared(dat.x[0][t]-dat.x[0][t-1],dat.x[1][t]-dat.x[1][t-1]);
  threshold=percentile(dx_squared,0.99);
  t01.push_back(0);
  for(int t=1;t<dat.T;t++){
    if(t!=1 && dx_squared[t-1]<threshold && (t==dat.T-1 || dx_squared[t]>=threshold))
      t01.push_back(t),n++;
    if((t==1 || dx_squared[t-1]>=threshold) && t!=dat.T-1 &&dx_squared[t]<threshold)
      t10.push_back(t);
    if(dx_squared[t]>=threshold)
      N++;
  }
  t10.push_back(dat.T);
}

void changepoints::write(char* filename){
  ofstream output(filename,ios::out | ios::app);
  output<<"N\t"<<N<<endl;
  output<<"n\t"<<n<<endl;
  output<<"t01\t";
  for(int i=0;i<n;i++)
    output<<t01[i]<<"\t";
  output<<endl;
  output<<"t10\t";
  for(int i=0;i<n;i++)
    output<<t10[i]<<"\t";
  output<<endl;
  output.close();
}

void changepoints::write(){
  cout<<"N\t"<<N<<endl;
  cout<<"n\t"<<n<<endl;
  cout<<"t01\t";
  for(int i=0;i<n;i++)
    cout<<t01[i]<<"\t";
  cout<<endl;
  cout<<"t10\t";
  for(int i=0;i<n;i++)
    cout<<t10[i]<<"\t";
  cout<<endl;
}

bool MCMC_chain::accept(double deltaL){
  exponential_distribution<double> dist(1.0);
  return (-dist(engine))<deltaL;
}


inline double data::get_dz_squared(const int t1,const int t2){
  if (t1==0)
    return normsquared(z[0][t2-1],z[1][t2-1]);
  return normsquared(z[0][t2-1]-z[0][t1-1],z[1][t2-1]-z[1][t1-1]);
}

inline double calc_loglik_up(int t1,int t2, params& theta, data& dat, integraltable& table){
  double dzsquared=dat.get_dz_squared(t1,t2);
  int dt=t2-t1;
  double a=0.5*pow(theta.sigmaz,4)/dzsquared*(1.0/pow(theta.sigma1,2)+dt/pow(theta.sigmaz,2));
  return theta.const_term_up-theta.prefactor*dt-0.5*(theta.d1+1.0)*log(dzsquared)+table(a,theta.d1);
}

inline double calc_loglik_down(int t1,int t2, params& theta, data& dat, integraltable& table){
  double dzsquared=dat.get_dz_squared(t1,t2);
  int dt=t2-t1;
  double a=0.5*pow(theta.sigmaz,4)/dzsquared*(1.0/pow(theta.sigma0,2)+dt/pow(theta.sigmaz,2));
  return theta.const_term_down-theta.prefactor*dt-0.5*(theta.d0+1.0)*log(dzsquared)+table(a,theta.d0);
}

inline double calc_loglik(data& dat, changepoints& C, params& theta, integraltable& table){
  double L=-0.5*dat.sum_dz_squared/(theta.sigmaz*theta.sigmaz);
  for(int k=0;k<C.n;k++)
    L+=calc_loglik_up(C.t01[k],C.t10[k],theta,dat,table);
  for(int k=0;k<C.n-1;k++)
    L+=calc_loglik_down(C.t10[k],C.t01[k+1],theta,dat,table);
  return L;
}

inline double calc_logprior(const params& theta, changepoints& C, const data& dat){
  int n=C.n;
  double lambda0=theta.lambda0, lambda1=theta.lambda1;
  double logpr=n*lambda0+(n-1)*lambda1+2.0*n*log(1.0-exp(-lambda0))+2.0*(n-1)*log(1.0-exp(-lambda1));
  for(int k=0;k<C.n;k++)
    logpr+=log(C.t10[k]-C.t01[k])-lambda1*(C.t10[k]-C.t01[k]);
  for(int k=0;k<C.n-1;k++)
    logpr+=log(C.t01[k+1]-C.t10[k])-lambda0*(C.t01[k+1]-C.t10[k]);
  return logpr;
}

void data::write_z(char* filename){
  ofstream output(filename,ios::out | ios::app);
  output<<T<<endl;
  for(int i=0;i<T;i++)
  	output<<z[0][i]<<"\t"<<z[1][i]<<endl;
  output.close();
}

inline double MCMC_chain::get_logprior(){
  return calc_logprior(theta, C, *dat);
}



inline double MCMC_chain::loglik_up(int t1, int t2){
  return calc_loglik_up(t1,t2,theta,*dat,*table);
}

inline double MCMC_chain::loglik_down(int t1, int t2){
  return calc_loglik_down(t1,t2,theta,*dat,*table);
}

inline double MCMC_chain::get_loglik(){
  return calc_loglik(*dat,C,theta,*table);
}

inline double MCMC_chain::get_logpost(){
  return get_loglik()+get_logprior();
}


bool MCMC_chain::make_step(const params& theta){
  bool smth_changed=true;
  int k=0;
  int pos=0;
  double comp=0.0;
  double deltaloglik=0.0;
  double deltalogprior=0.0;
  double logprior_new=logprior;
  int step_type=get_step_type(engine);
  int n=C.n;
  const double lambda0=theta.lambda0, lambda1=theta.lambda1;
  const double log_norm=lambda0+lambda1+2.0*log(1.0-exp(-lambda0))+2.0*log(1.0-exp(-lambda1));
  switch(step_type){
  case 0: //select a t01, make it t01-1
  if (n>=2){
    k=uniform_int_distribution<int>{1,n-1}(engine);
    pos=C.t01[k]-1;
    if (C.t01[k]-C.t10[k-1]==1){ //if this eliminates a changepoint
      deltaloglik=loglik_up(C.t01[k-1],C.t10[k])-loglik_up(C.t01[k-1],C.t10[k-1])-loglik_down(C.t10[k-1],C.t01[k])-loglik_up(C.t01[k],C.t10[k]);
      comp=log(C.t10[k]-C.t01[k-1]-2);
      deltalogprior=-log_norm+lambda0-lambda1+log(C.t10[k]-C.t01[k-1])-log(C.t10[k-1]-C.t01[k-1])-log(C.t10[k]-C.t01[k]);
    }
    else {
      deltaloglik=loglik_up(pos,C.t10[k])+loglik_down(C.t10[k-1],pos)-loglik_up(C.t01[k],C.t10[k])-loglik_down(C.t10[k-1],C.t01[k]);
      deltalogprior=lambda0-lambda1+log(C.t10[k]-pos)+log(pos-C.t10[k-1])-log(C.t10[k]-C.t01[k])-log(C.t01[k]-C.t10[k-1]);
    }
  }else smth_changed=false;
  break;
  case 1: // select a C.t01, make it C.t01+1
  if (n>=2){
    k=uniform_int_distribution<int>{1,n-1}(engine);
    pos=C.t01[k];
    if (C.t10[k]-C.t01[k]==1){
      if (k!=n-1){
        deltaloglik=loglik_down(C.t10[k-1],C.t01[k+1])-loglik_down(C.t10[k-1],C.t01[k])-loglik_up(C.t01[k],C.t10[k])-loglik_down(C.t10[k],C.t01[k+1]);
        comp=-log(n-1)+log(n-2)+log(C.t01[k+1]-C.t10[k-1]-2);
        deltalogprior=-log_norm+lambda1-lambda0+log(C.t01[k+1]-C.t10[k-1])-log(C.t01[k]-C.t10[k-1])-log(C.t01[k+1]-C.t10[k]);
      }
      else smth_changed=false;
    }
    else{
      deltaloglik=loglik_up(pos+1,C.t10[k])+loglik_down(C.t10[k-1],pos+1)-loglik_up(C.t01[k],C.t10[k])-loglik_down(C.t10[k-1],C.t01[k]);
      deltalogprior=lambda1-lambda0+log(C.t10[k]-pos-1)+log(pos+1-C.t10[k-1])-log(C.t10[k]-C.t01[k])-log(C.t01[k]-C.t10[k-1]);
    }
  }else smth_changed=false;
  break;
  case 2: //select a C.t10 and make it C.t10+1
    if (n>=2){
      k=uniform_int_distribution<int>{0,n-2}(engine);
      pos=C.t10[k];
      if(C.t01[k+1]-C.t10[k]==1){ //if it eliminates a changepoint
        deltaloglik=loglik_up(C.t01[k],C.t10[k+1])-loglik_up(C.t01[k],C.t10[k])-loglik_down(C.t10[k],C.t01[k+1])-loglik_up(C.t01[k+1],C.t10[k+1]);
        comp=log(C.t10[k+1]-C.t01[k]-2);
        deltalogprior=-log_norm-lambda1+lambda0+log(C.t10[k+1]-C.t01[k])-log(C.t10[k]-C.t01[k])-log(C.t10[k+1]-C.t01[k+1]);
      }else{
        deltaloglik=loglik_up(C.t01[k],pos+1)+loglik_down(pos+1,C.t01[k+1])-loglik_up(C.t01[k],C.t10[k])-loglik_down(C.t10[k],C.t01[k+1]);
        deltalogprior=lambda0-lambda1+log(pos+1-C.t01[k])+log(C.t01[k+1]-pos-1)-log(C.t10[k]-C.t01[k])-log(C.t01[k+1]-C.t10[k]);
      }
    }else smth_changed=false;
  break;
  case 3: //  select a C.t10 and make it C.t10-1
    if (n>=2){
      k=uniform_int_distribution<int>{0,n-2}(engine);
      pos=C.t10[k]-1;
      if (C.t10[k]-C.t01[k]==1){ //if it eliminates a changepoint
        if (k!=0){
          deltaloglik=loglik_down(C.t10[k-1],C.t01[k+1])-loglik_down(C.t10[k-1],C.t01[k])-loglik_up(C.t01[k], C.t10[k])-loglik_down(C.t10[k],C.t01[k+1]);
          comp=-log(n-1)+log(n-2)+log(C.t01[k+1]-C.t10[k-1]-2);
          deltalogprior=lambda1-lambda0-log_norm+log(C.t01[k+1]-C.t10[k-1])-log(C.t01[k]-C.t10[k-1])-log(C.t01[k+1]-C.t10[k]);
        }else smth_changed=false;
      }else{
        deltaloglik=loglik_up(C.t01[k], pos)+loglik_down(pos,C.t01[k+1])-loglik_up(C.t01[k],C.t10[k])-loglik_down(C.t10[k], C.t01[k+1]);
        deltalogprior=lambda1-lambda0+log(pos-C.t01[k])+log(C.t01[k+1]-pos)-log(C.t10[k]-C.t01[k])-log(C.t01[k+1]-C.t10[k]);
      }
    }else smth_changed=false;
    break;
  case 4: //add a new pair by adding a 0 in a seq of 1s
    k=uniform_int_distribution<int>{0,n-1}(engine);
    if(C.t10[k]-C.t01[k]>2){
      pos=uniform_int_distribution<int>{C.t01[k]+1,C.t10[k]-2}(engine);
      deltaloglik=loglik_up(C.t01[k],pos)+loglik_down(pos,pos+1)+loglik_up(pos+1,C.t10[k])-loglik_up(C.t01[k],C.t10[k]);
      comp=-log(C.t10[k]-C.t01[k]-2);
      deltalogprior=log_norm+lambda1-lambda0+log(pos-C.t01[k])+log(C.t10[k]-pos-1)-log(C.t10[k]-C.t01[k]);
    } else smth_changed=false;
    break;
  case 5: //add a new pair by adding a 1 in a seq of 0s
    if (n>=2){
      k=uniform_int_distribution<int>{0,n-2}(engine);// changed 0 to 1
      if(C.t01[k+1]-C.t10[k]>2){
        pos=uniform_int_distribution<int>{C.t10[k]+1,C.t01[k+1]-2}(engine);
        deltaloglik=loglik_down(C.t10[k],pos)+loglik_up(pos,pos+1)+loglik_down(pos+1,C.t01[k+1])-loglik_down(C.t10[k],C.t01[k+1]);
        comp=-log(n-1)-log(C.t01[k+1]-C.t10[k]-2)+log(n);
        deltalogprior=log_norm+lambda0-lambda1+log(pos-C.t10[k])+log(C.t01[k+1]-pos-1)-log(C.t01[k+1]-C.t10[k]);
      } else smth_changed=false;
    } else smth_changed=false;
    break;
  }
  if (smth_changed && accept(deltalogprior+beta*deltaloglik-comp)){
    C.execute_step(step_type,k,pos);
    logprior+=deltalogprior;
    loglik+=deltaloglik;
    logposterior=logprior+loglik;
    if(logposterior>logpostbest)
      logpostbest=logposterior;
    nsteps++;
    return true;
  }
  return false;
}

void changepoints::execute_step(int step_type,int k, int pos){
  switch(step_type){
  case 0: //select a t01, make it t01-1
  if (t01[k]-t10[k-1]==1){  //if it eliminates a changepoint pair
  	t01.erase(t01.begin()+k);
  	t10.erase(t10.begin()+k-1);
  	n--;
  	N++;
  }else{
  	t01[k]--;
  	N++;
  }
  break;
  case 1://select a t01, make it a t01+1
  if (t10[k]-t01[k]==1){ //if this eliminates a changepoint
    t01.erase(t01.begin()+k);
    t10.erase(t10.begin()+k);
    n--;
    N--;
  }else{
    t01[k]++;
    N--;
  }
  break;
  case 2://select a t10 and make it t10+1
  if (t01[k+1]-t10[k]==1){ //if it eliminates a changepoint
  	t01.erase(t01.begin()+k+1);
  	t10.erase(t10.begin()+k);
    n--;
    N++;
  }else{
  	t10[k]++;
  	N++;
  }
  break;
  case 3://  select a t10 and make it t10-1
  if (t10[k]-t01[k]==1){ //if it eliminates a changepoint
  	t01.erase(t01.begin()+k);
  	t10.erase(t10.begin()+k);
  	n--;
  	N--;
  }else{
  	t10[k]--;
  	N--;
  }
  break;
  case 4: //add a new pair by adding a 0 in a seq of 1s
    t01.insert(t01.begin()+k+1,pos+1);
    t10.insert(t10.begin()+k,pos);
    n++;
    N--;
    break;
  case 5: //add a new pair by adding a 1 in a seq of 0s
    t01.insert(t01.begin()+k+1,pos);
    t10.insert(t10.begin()+k+1,pos+1);
    n++;
    N++;
    break;
  }
}

bool MCMC_chain::swap_chains(MCMC_chain* other){
  swaps_tried++;
  if(accept((other->beta-beta)*(loglik-other->loglik))){
    swaps_acc++;
    swap(beta,other->beta);
    swap(swaps_acc,other->swaps_acc);
    swap(swaps_tried,other->swaps_tried);
    return true;
  }
  return false;
}

void MCMC_chain::change_params(const params& theta_new){
  theta=theta_new;
  logprior=get_logprior();
  loglik=get_loglik();
  logposterior=loglik+logprior;
}

void parallel_tempered_chain::update_params(){
  for(int i=0;i<Nchains;i++)
    if(chain[i])
      chain[i]->change_params(theta);
}

void parallel_tempered_chain::sweep(){
  int chi;
  for(int i=0;i<dat.T;i++){
    if(Nchains>1 && i%10==0){
      chi=random_chain(engine);
      if(chain[chi]->swap_chains(chain[chi+1]))
        swap(chain[chi],chain[chi+1]);
    }
    else for(int j=0;j<Nchains;j++)
      chain[j]->make_step(theta);
  }
}

void parallel_tempered_chain::write_posteriors(int n,char* filename_posts,char* filename_swaps){
  ofstream output_swaps(filename_swaps,ios::out);
  ofstream output_posts(filename_posts,ios::out);
  for(int i=0;i<n;i++){
  	for(int j=0;j<Nchains;j++)
      chain[j]->swaps_acc=chain[j]->swaps_tried=0;
  	sweep();
  	for(int j=0;j<Nchains;j++)
      output_swaps<<chain[j]->swaps_acc<<"\t"<<chain[j]->swaps_tried<<"\t";
    output_swaps<<endl;
  	for(int j=0;j<Nchains;j++)
  	  output_posts<<chain[j]->logposterior<<"\t";
  	output_posts<<endl;
  }
  output_swaps.close();
  output_posts.close();
}

void parallel_tempered_chain::sweep(int nsweeps){
  for(int i=0;i<nsweeps;i++){
    sweep();
    cout<<"sweep\t"<<i<<"/"<<nsweeps<<":\t"<<chain[neo]->logposterior<<endl;
  }
}

void parallel_tempered_chain::initialize_changepoints(){
  C_samples.emplace_back(dat);
  nsamples=1;
}

void parallel_tempered_chain::sample_changepoints(int n){
  C_samples.clear();
  nsamples=n;
  for(int i=0;i<nsamples;i++){
    sweep();
    C_samples.push_back(chain[neo]->C);
    cout<<"sample\t"<<i<<"/"<<nsamples<<":\t"<<chain[neo]->logposterior<<endl;
  }
}

void parallel_tempered_chain::write_samples(char* filename){
  ofstream output;
  for(int i=0;i<nsamples;i++){
    C_samples[i].write(filename);
    output.open(filename,ios::out | ios::app);
    output<<"logpost\t"<<get_logpost(theta,C_samples[i])<<endl;
    //output<<"logpost\t"<<get_logpost(theta,C_samples[i])<<"logpostreal\t"<<get_logpost(theta,C_real)<<endl;
    output.close();
  }
}

void parallel_tempered_chain::estimate_dsigma_up(){
  double d,sigma,d_best,sigma_best;
  params theta_new=theta;
  double L,L_best;
  vector<double> sigma_est,d_est;
  cout<<"estimate dsigma_up"<<endl;
  for(int i=0;i<nsamples;i++){
    for(int di=0;di<100;di++){
  	  d=0.1+di*(5.0-0.1)/100.0;
  	  for(int sgi=0;sgi<100;sgi++){
        sigma=pow(10.0, -4.0+4.0/100.0*sgi);
  	    theta_new.set_d_sigma_up(d,sigma);
        L=get_logpost(theta_new,C_samples[i]);
        if((di==0 && sgi==0) || L>L_best)
          sigma_best=sigma,d_best=d,L_best=L;
  	  }
  	}
  	sigma_est.push_back(sigma_best);
  	d_est.push_back(d_best);
  }
  sort(sigma_est.begin(),sigma_est.end());
  sort(d_est.begin(),d_est.end());
  theta.set_d_sigma_up(d_est[nsamples/2],sigma_est[nsamples/2]);
  update_params();
  cout<<"d1="<<theta.d1<<"\tsigma1="<<theta.sigma1<<endl;
}

void parallel_tempered_chain::estimate_sigma_down(){
  double sigma,sigma_best;
  params theta_new=theta;
  double L,L_best;
  vector<double> sigma_est;
  cout<<"estimate sigma_down"<<endl;
  for(int i=0;i<nsamples;i++){
    for(int sgi=0;sgi<100;sgi++){
  	  sigma=pow(10.0, -4.0+4.0/100.0*sgi);
  	  theta_new.set_sigma_down(sigma);
      L=get_logpost(theta_new,C_samples[i]);
      //cout<<"L="<<L<<endl;
      if(sgi==0 || L>L_best)
        sigma_best=sigma,L_best=L;
    }
    sigma_est.push_back(sigma_best);
  }
  sort(sigma_est.begin(),sigma_est.end());
  theta.set_sigma_down(sigma_est[nsamples/2]);
  update_params();
  cout<<"d0="<<theta.d0<<"\tsigma0="<<theta.sigma0<<endl;
}

void data::get_residual(changepoints& C, double* x_res[2]){
  double a;
  double pred;
  int i,j;
  for(int dim=0;dim<2;dim++){
    i=0,j=0;
    pred=x[dim][0];
    for(int t=0;t<T;t++){
      if(C.t01[i]==t)
        i++,a=(x[dim][C.t10[j]]-x[dim][t])/(C.t10[j]-t);
      if(C.t10[j]==t)
        j++,a=(x[dim][C.t01[i]]-x[dim][t])/(C.t01[i]-t);
      x_res[dim][t]=x[dim][t]-pred;
      pred+=a;
    }
  }
}

void parallel_tempered_chain::estimate_sigmazx(){
  const int nmax=20;
  double* x_res[2];
  double Corr;
  const double mean_n=nmax/2.0;
  const double mean_nn=nmax*(2*nmax-1)/6.0;
  double mean_Corr;
  double mean_Corrn;
  double slope,intercept;
  vector<double> sigmaz_est,sigmax_est;
  cout<<"estimate sigma_zx"<<endl;
  x_res[0]=new double[dat.T];
  x_res[1]=new double[dat.T];
  for(int i=0;i<nsamples;i++){
    dat.get_residual(C_samples[i],x_res);
    mean_Corr=0.0;
    mean_Corrn=0.0;
    for(int n=1;n<nmax;n++){
      Corr=0.0;
      for(int t=n;t<dat.T;t++)
        Corr+=normsquared(x_res[0][t]-x_res[0][t-n],x_res[1][t]-x_res[1][t-n]);
      Corr/=(dat.T-n);
      mean_Corr+=Corr/(nmax-1);
      mean_Corrn+=n*Corr/(nmax-1);
    }
    slope=(mean_Corrn-mean_Corr*mean_n)/(mean_nn-mean_n*mean_n);
    intercept=(mean_Corr-slope*mean_n);
    if (abs(intercept-0.0)<0.00001)
       intercept=0.00001;
    sigmaz_est.push_back(sqrt(slope/2.0));
    sigmax_est.push_back(sqrt(max(0.00001,intercept)/4.0));
  }
  sort(sigmax_est.begin(),sigmax_est.end());
  sort(sigmaz_est.begin(),sigmaz_est.end());
  theta.set_sigmaz(sigmaz_est[nsamples/2]);
  theta.set_sigmax(sigmax_est[nsamples/2]);
  update_params();
  cout<<intercept<<"\tsigma_z="<<theta.sigmaz<<"\tsigma_x="<<theta.sigmax<<endl;
  delete x_res[0];
  delete x_res[1];
}

void parallel_tempered_chain::estimate_z(){
  cout<<"estimate z"<<endl;
  dat.kalman_filter(theta);
  update_params();
}

double parallel_tempered_chain::get_logprior(params& p, changepoints& C){
  return calc_logprior(p,C,dat);
}

double parallel_tempered_chain::get_loglik(params& p, changepoints& C){
  return calc_loglik(dat,C,p,table);
}

double parallel_tempered_chain::get_logpost(params& p, changepoints& C){
  return get_loglik(p,C)+get_logprior(p,C);
}

void parallel_tempered_chain::estimate_allthethings(char* filename_params, char* filename_samples){
  for(int i=0;i<6;i++){
    theta.write(filename_params);
    //cout<<"logpostreal\t"<<get_logpost(theta,C_real)<<"\tlogpostinitial\t"<<chain[neo]->get_logpost()<<endl;
    cout<<"logpostreal\t"<<"\tlogpostinitial\t"<<chain[neo]->get_logpost()<<endl;
    sweep(40);
    sample_changepoints(40);
    write_samples(filename_samples);
    estimate_dsigma_up();
    estimate_sigma_down();
    estimate_sigmazx();
    estimate_z();
  }
}

parallel_tempered_chain::~parallel_tempered_chain(){
  for(int i=0;i<Nchains;i++)
    delete chain[i];
}

parallel_tempered_chain::parallel_tempered_chain(double beta[], char* filename_dat, char* filename_table):
random_chain(0,Nchains-2), dat(filename_dat), table(filename_table){
  for(int i=0;i<Nchains;i++)
    chain[i]=NULL;
  //unsigned int seed=unsigned(time(0));
  unsigned int seed=1437768682;
  engine.seed(seed);
  cout<<"seed:\t"<<seed<<endl;
  /*
  //To start with the real C
  C_samples.push_back(C_real);
  nsamples=1;
  */
  //To start from the data
  initialize_changepoints();
  estimate_sigmazx();
  estimate_z();
  for(int i=0;i<Nchains;i++){
    chain[i]=new MCMC_chain(beta[i],&dat,theta,C_samples[0],&table);
    chain[i]->engine.seed(engine());
  }
}


int main(int argc, char* argv[]){
  double beta[82]={0.0,0.0005,0.001,0.0015,0.002,0.003, 0.004,0.005,0.006,0.007,0.0085,0.01,0.0115,
                   0.013,0.015,0.017,0.0185,0.02,0.023, 0.025,0.027,0.03,0.033,0.037,
                   0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.09,0.1,0.11,0.12,
                   0.13,0.14,0.15,0.17,0.185,0.2,0.23, 0.25,0.27,0.3,0.35,0.4,0.45,0.5,
                   0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.5, 1.6,1.7,1.8,1.9, 2.0, 2.1, 2.2,
                52.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.25, 4.5};
  if(argc<6){
    cout<<"I need at least 6 arguments:"<<endl
    <<"x_input_file\tintegral_table\tparams_out\tchangepoints_out\tparal_temp\tchain_swaps\t"<<endl;
    return 1;
  }
  parallel_tempered_chain mychains(beta, argv[1], argv[2]);
  /*
  argv[1]: filename for data
  argv[2]: filename for table
  */
  mychains.estimate_allthethings(argv[3],argv[4]);
  //if(argc==7)
  mychains.write_posteriors(82,argv[5],argv[6]);
  return 0;
}
