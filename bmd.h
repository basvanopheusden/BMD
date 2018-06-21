#ifndef CHANGEPOINTS_H_INCLUDED
#define CHANGEPOINTS_H_INCLUDED

#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;

struct integraltable{
  double* table;
  double* a_upper;
  double* a_lower;
  int Nd;
  int Na;
  double d_step_size;
  integraltable(char*);
  ~integraltable();
  int get_a_ind(const double);
  int get_d_ind(const double);
  double get_a(const int);
  double get_d(const int);
  double operator()(const double,const double);
};

struct params{
    double  lambda0=0.004, lambda1=0.1; // parameters of the prior distributions (gamma) over durations
  long double d0=1.0;
  long double d1,sigma0,sigma1,sigmaz,sigmax;
  double const_term_up;
  double const_term_down;
  double prefactor;
  mt19937_64 engine;
  params();
  params(char*);
  params(const double, const double, const double, const double, const double);
  void write(char*);
  void write();
  void set_d_sigma_up(const double, const  double);
  void set_sigma_down(const double);
  void set_sigmaz(const double);
  void set_sigmax(const double);
};

struct changepoints;

struct data{
  double* x[2];
  double* z[2];
  double* zf[2];
  int T;
  double sum_dz_squared;
  data(char*);
  ~data();
  void write_z(char*);
  void kalman_filter(const params& );
  double get_dz_squared(int,int);
  void get_residual(changepoints&, double*[2]);
};

struct changepoints{
  vector<int> t01;
  vector<int> t10;
  int n,N;
  vector<bool> get_vec();
  changepoints(int);
  changepoints(data&);
  changepoints(char*,int);
  void write(char*);
  void write();
  void execute_step(int,int,int); //makes a step with a particular type, changepoint index k and position
};

struct MCMC_chain{
  double beta;
  data* dat;
  params theta;
  changepoints C;
  integraltable* table;
  double logposterior;
  double logprior;
  double loglik;
  double logpostbest;
  int nsteps;
  int swaps_acc;
  int swaps_tried;
  mt19937_64 engine;
  discrete_distribution<int> get_step_type;
  MCMC_chain(double,data*,params&,changepoints&,integraltable*);
  bool make_step(const params&); //generates a step type and k, accepts with Metropolis probability
    //and returns true if successful
  bool accept(double);//implements the metropolis acceptance probability
  bool swap_chains(MCMC_chain*);
  double get_logprior();
  double get_logprior(const int,const int);
  double loglik_up(const int,const int);
  double loglik_down(const int,const int);
  double get_loglik();
  double get_logpost();
  void change_params(const params&);
};

class parallel_tempered_chain{
  public:
    //static const int Nchains=82;
    static const int Nchains=1;
    MCMC_chain* chain[Nchains];
    uniform_int_distribution<int> random_chain;
    //static const int neo=55;
    static const int neo=0;
    mt19937_64 engine;
    data dat;
    params theta;
    //changepoints C_real; //if we want the Creal input file
    params theta_real;
    integraltable table;
    double logpostreal;
    vector<changepoints> C_samples;
    int nsamples;
    //parallel_tempered_chain(double*,char*,char*,char*); //if we want the Creal input file
    parallel_tempered_chain(double*,char*,char*);
    ~parallel_tempered_chain();
    double get_loglik(params&,changepoints&);
    double get_logprior(params&,changepoints&);
    double get_logpost(params&,changepoints&);
    void update_params();
    void sweep();
    void sweep(int);
    void estimate_allthethings(char*,char*);
    void sample_changepoints(int);
    void initialize_changepoints();
    void estimate_z();
    void estimate_sigmazx();
    void estimate_dsigma_up();
    void estimate_sigma_down();
    void write_posteriors(int,char*, char*);
    void write_changepoints(char*);
    void write_samples(char*);
};

#endif // CHANGEPOINTS_H_INCLUDED
