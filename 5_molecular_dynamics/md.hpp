#ifndef MD_HPP
#define MD_HPP
#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <complex>
#include <bitset>
#include <time.h>
using std::ofstream;
#include <random>
#include <omp.h>


random_device rd{};
mt19937 gen{rd()};

// \tau mora biti večkratnik dt!!
struct mch{ // molecular chain
  int   N; //number of chain links
  double dt,
         TL,
         TR,
         lambda,
         tdead,
         tau,
         tmax,
         tcur;


  VectorXd state; // [q0,p0,q1,p1,...,qN-1,pN-1,KSI_L,KSI_R]
  VectorXd Tavg;
  VectorXd Javg;
  const VectorXd mass;
  //VectorXd stateb; // state_before (one time step before)

  mch(int N_,double dt_, double TL_, double TR_, double lambda_, double tdead_,
  double tmax_,double tau_, VectorXd *ptr_mass=NULL, VectorXd *ptr_state=NULL ):
  N(N_),dt(dt_), TL(TL_),TR(TR_), lambda(lambda_), tdead(tdead_), tau(tau_),
  tmax(tmax_), mass((ptr_mass != NULL) ? *ptr_mass : VectorXd::Ones(N)){

  if (ptr_state!=NULL){
     state = *ptr_state;
  }
  else {
    srand (12);
    state = MatrixXd::Random(2*N+2,1).cwiseAbs();


    }
  //stateb = VectorXd::Zero(2*N+2,1);
  Tavg = VectorXd::Zero(N);
  Javg = VectorXd::Zero(N);
  tcur = 0;
  }
};

// NOSE-HOOVER
VectorXd Func(VectorXd  statt, mch &ch){
  VectorXd newstate(2*ch.N +2 );
  //ksiL:
  newstate(2*ch.N) = (1/ch.tau) * ((pow(statt(1), 2) / ch.mass(0)) - ch.TL);

  //ksiR:
  newstate(2*ch.N+1) = (1/ch.tau) * ((pow(statt(2*ch.N-1), 2) / ch.mass(ch.N-1)) - ch.TR);
  //robni enacbi:
  //levi rob:

  newstate(0) = statt(1)/ch.mass(0);
  newstate(1) = -(2*statt(0) + 4*ch.lambda*pow(statt(0), 3) -
  statt(2)) - statt(1)*statt(2*ch.N);
  //desni:
  newstate(2*ch.N-2) = statt(2*ch.N-1)/ch.mass(ch.N-1);
  newstate(2*ch.N-1) = -(2*statt(2*ch.N-2) + 4*ch.lambda*pow(statt(2*ch.N-2), 3)
  - statt(2*ch.N-4)) - statt(2*ch.N + 1)*statt(2*ch.N -1);
  // ostali:

  #pragma omp parallel
  {
  #pragma omp for
  for (size_t i = 1; i < (ch.N-1); i++) {
    newstate(2*i) = statt(2*i+1) / ch.mass(i);
    newstate(2*i+1) = -(3*statt(2*i) + 4*ch.lambda*pow(statt(2*i), 3) -
    statt(2*(i-1))  - statt(2*(i+1)));

    }
  }
  std::cout << newstate << '\n';
  return newstate;
}

void RK4_one_step(mch &ch){
  VectorXd k1, k2, k3, k4;

  k1 = ch.dt *Func(ch.state, ch);
  k2 = ch.dt*Func(ch.state + 0.5*k1, ch);
  k3 = ch.dt*Func(ch.state + 0.5*k2, ch);
  k4 = ch.dt*Func(ch.state + k3, ch);
  ch.state += (1.0/6.0)*k1 + (1.0/3)*k2 +  (1.0/3)*k3 +(1.0/6.0)*k4;
}

void relaxation(mch &ch, std::string metoda){
  if (metoda == "hoover"){
    unsigned int N = (unsigned int) ch.tdead/ch.dt;
    for (size_t i = 0; i < N; i++)
    {
      RK4_one_step(ch);
    }

  }
  else if (metoda=="maxwell"){
    srand(time(NULL));
    normal_distribution<> p0{0, sqrt(ch.mass(0)*ch.TL)};
    normal_distribution<> pn{0, sqrt(ch.mass(ch.N-1)*ch.TR)};

    auto N = (int) ch.tdead/ch.dt;
    double strobotime = ch.tau;
    double tcurd = 0;
    for (size_t i = 0; i < N; i++) {
      if (abs(tcurd - strobotime) < pow(10, -7)){
        ch.state(1) = p0(gen);
        ch.state(2*ch.N-1) = pn(gen);
        strobotime += ch.tau;
        }
      else{
        RK4_one_step(ch);
    }
      tcurd += ch.dt;
    }
  }
  else {
    std::cout << "izbrana metoda ne obstaja" << '\n';
    exit(1);
  }
}

VectorXd calculateJ(mch &ch){
  VectorXd Jvec(ch.N);
  for (size_t i = 1; i <(ch.N-1); i++) {
    Jvec(i) = -0.5*(ch.state(2*(i+1)) - ch.state(2*(i-1)))*ch.state(2*i +1 );
  }
  Jvec(0) = Jvec(1);
  Jvec(ch.N-1) = Jvec(ch.N-2);
  return Jvec;
}


// using trapez formula to average:
void add_to_avg(mch &ch){

  if (ch.tcur==0 || abs(ch.tcur-(ch.tmax-ch.dt)) < pow(10, -10)) {
    for (size_t i = 0; i < ch.N; i++) {
      ch.Tavg(i) += pow(ch.state(2*i+1),2);
    }
    ch.Javg += calculateJ(ch);

  }
  else{
    for (size_t i = 0; i < ch.N; i++) {
      ch.Tavg(i) += 2.0*pow(ch.state(2*i+1),2);
    }
    ch.Javg += 2*calculateJ(ch);
  }
}

void RK4_iteration(mch &ch){
  relaxation(ch, "hoover");

  auto N = (int) ch.tmax/ch.dt;
  add_to_avg(ch);
  for (size_t i = 0; i < N-1; i++) {
    RK4_one_step(ch);
    ch.tcur += ch.dt;
    add_to_avg(ch);
  }
  ch.Tavg = ch.dt/(2*ch.tcur)* ch.Tavg;
  ch.Javg = ch.dt/(2*ch.tcur)* ch.Javg;
}


// MAXWELL BATH:
VectorXd Func_mawxell(VectorXd  statt, mch &ch){
  VectorXd newstate(2*ch.N +2 );
  //robni enacbi:
  //levi rob:
  newstate(0) = statt(1)/ch.mass(0);
  newstate(1) = -(2*statt(0) + 4*ch.lambda*pow(statt(0), 3) -
  statt(2));
  //desni:
  newstate(2*ch.N-2) = statt(2*ch.N-1)/ch.mass(ch.N-1);
  newstate(2*ch.N-1) = -(2*statt(2*ch.N-2) + 4*ch.lambda*pow(statt(2*ch.N-2), 3)
  - statt(2*ch.N-4));
  // ostali:
  for (size_t i = 1; i < (ch.N-1); i++) {
    newstate(2*i) = statt(2*i+1) / ch.mass(i);
    newstate(2*i+1) = -(3*statt(2*i) + 4*ch.lambda*pow(statt(2*i), 3) -
    statt(2*(i-1))  - statt(2*(i+1)));
  }

  return newstate;
}

void RK4_one_step_max(mch &ch){
  VectorXd k1, k2, k3, k4;
  k1 = ch.dt *Func_mawxell(ch.state, ch);
  k2 = ch.dt*Func_mawxell(ch.state + 0.5*k1, ch);
  k3 = ch.dt*Func_mawxell(ch.state + 0.5*k2, ch);
  k4 = ch.dt*Func_mawxell(ch.state + k3, ch);
  ch.state += (1.0/6.0)*k1 + (1.0/3)*k2 +  (1.0/3)*k3 +(1.0/6.0)*k4;
}


void RK4_maxwell_prop(mch &ch){
  relaxation(ch, "maxwell");

  normal_distribution<> p00{0, sqrt(ch.mass(0)*ch.TL)};
  normal_distribution<> pnn{0, sqrt(ch.mass(ch.N-1)*ch.TR)};

  auto N = (int) ch.tmax/ch.dt;
  add_to_avg(ch);
  double strobotime=ch.tau;
  for (size_t i = 0; i < N-1; i++) {
    ch.tcur += ch.dt;
    if (abs(ch.tcur - strobotime) < pow(10, -7)){
      ch.state(1) = p00(gen);
      ch.state(2*ch.N-1) = pnn(gen);
      strobotime += ch.tau;
      add_to_avg(ch);
      }
    else{
      RK4_one_step(ch);
      add_to_avg(ch);
    }
  }
  ch.Tavg = ch.dt/(2*ch.tcur)* ch.Tavg;
  ch.Javg = ch.dt/(2*ch.tcur)* ch.Javg;

}


















#endif
