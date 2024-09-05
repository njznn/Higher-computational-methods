#ifndef QMC_HPP
#define QMC_HPP
#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <complex>
#include <time.h>
using std::ofstream;
#include <random>
#include <omp.h>

random_device rd{};
mt19937 gen{rd()};


struct qmc {
  int    M, //number of imag time slices
         N; // number of chain links
  double B, // beta
         eps, // model parameter for accepting qmc moves;
         lam, // quartic potential parameter
         Vavg, // avg potential energy;
         Tavg, // avg kinetic energy;
         Eavg,
         delez;


  int nrelax; // number of mc relaxation steps
  int nsample; // -||- sampling steps
  int Nacept = 0; // number of accepting moves//

  MatrixXd state; // rows represent chain (N), columns time slices (M)

  qmc(double M_, double B_, double N_, double eps_, double lam_, int nrelax_
  , int nsample_, double seed_ , std::string potencial): M(M_), B(B_), N(N_),
  eps(eps_), lam(lam_), nrelax(nrelax_), nsample(nsample_){
    if (potencial == "kvad" && potencial != "kvar"){
      lam = 0;
    }
    else if (potencial != "kvar"){
      std::cout << "Izberi pravo obliko potenciala!" << '\n';
      exit(1);
    }
    gen.seed(seed_);
    srand(seed_);
    state =  MatrixXd::Random(M_,N_);
  }
};

// 0.5*q_^2 + \lambda * q_^4, j is time slice!
double Vslice(qmc & obj, int j){
  double pot = 0;
  for (size_t i = 0; i < obj.N; i++) {
    pot += 0.5*pow(obj.state(j, i), 2) + obj.lam * pow(obj.state(j, i), 4);

  }
  return pot;
}

//interaction between slices-kinetic energy:
//0.5*|q_(j+1)-q_(j)|^2
double Tslice(qmc & obj, int j){
  return 0.5*(obj.state((j+1)% obj.M, all)- obj.state(j%obj.M, all)).squaredNorm();
}

double Epoli(qmc & obj){
  double E = 0;
  for (size_t j = 0; j < obj.M; j++) {
    E += 1.0*(obj.M/obj.B)* Tslice(obj, j) + 1.0*(obj.B/obj.M)*Vslice(obj, j);
  }
  return E;
}

double Pjjnext(qmc & obj, int j){
  return(exp(-1.0*(obj.M/obj.B)* Tslice(obj, j) - 1.0*(obj.B/obj.M)*Vslice(obj, j) ));
}

double Eavgcalc(qmc &obj){
  double E = 1.0*obj.M/(2*obj.B);
  double A=0;
  double B=0;
  for (size_t i = 0; i < obj.M; i++) {
    A += Tslice(obj, i);
    B += Vslice(obj, i);
  }
  double res = (E - (1.0*obj.M/(pow(obj.B, 2) * obj.N))*A + (1.0/(obj.M*obj.N)) *B );

  return res;
}
double Vavgcalc(qmc &obj){
  double V = 0;
  for (size_t i = 0; i < obj.M; i++) {
    V += Vslice(obj, i);
  }
  return V/obj.M;
}

double Tavgcalc(qmc &obj){
  double E0 = 1.0*obj.M/(obj.B*2);
  for (size_t j = 0; j < obj.M; j++) {
    E0 -= obj.M/(pow(obj.B,2)) * Tslice(obj, j);
  }
  return E0;
}




void qmc_step(qmc &obj, int j, double K, VectorXd   ksirand){
  MatrixXd oldslice = obj.state(j, all);
  double Pj, Pjmin, Pj_, Pjmin_;
  Pj = Pjjnext(obj, j);
  Pjmin = (j==0) ? Pjjnext(obj, obj.M-1): Pjjnext(obj, j-1);
  if (obj.N>1){
    ksirand.normalize();
  }
  ksirand *= obj.eps;
  obj.state(j, all) += ksirand;
  Pj_ = Pjjnext(obj, j);
  Pjmin_ = (j==0) ? Pjjnext(obj, obj.M-1): Pjjnext(obj, j-1);
  //total probability
  double P = (Pj_*Pjmin_)/(Pj*Pjmin);
  if (P >= 1.0){
    obj.Nacept +=1;
  }
  else{
    if (K < P){
      obj.Nacept +=1;
    }
    else{
      obj.state(j, all) = oldslice;
    }

  }
}



void relax_and_sample(qmc &obj){
  std::normal_distribution<double> dist(0,1.0/obj.N);
  std::uniform_int_distribution<int> jpool(0,obj.M-1);
  std::uniform_real_distribution<double> weight(0.0,1.0);
  VectorXd hist = VectorXd::Zero(obj.M);
  VectorXd ksirand = VectorXd::Zero(obj.N);
  for (size_t j = 0; j < obj.nrelax; j++) {
    for (size_t k = 0; k < obj.N; k++) {
      ksirand(k) = dist(gen);
    }
    auto jaj = jpool(gen);
    qmc_step(obj, jaj, weight(gen), ksirand);
    hist +=obj.state;
  }

  int Navg = 0;
  double Vavg = 0;
  double Eavg = 0;
  double Tavg = 0;
  for (size_t i = 0; i < obj.nsample; i++) {
    for (size_t k = 0; k < obj.N; k++) {
      ksirand(k) = dist(gen);
    }
    qmc_step(obj, jpool(gen), weight(gen), ksirand);
    if (i%3000 == 0){
      Navg +=1;
      Eavg += Eavgcalc(obj);
      Vavg += Vavgcalc(obj);
      Tavg += Tavgcalc(obj);
    }
  }

  obj.Eavg = Eavg/Navg;
  obj.Vavg = Vavg/Navg;
  obj.Tavg = Tavg/Navg;
  obj.delez = (double) obj.Nacept/(obj.nsample + obj.nrelax);

}

void relax_and_hist(qmc &obj, std::string ime_dat){
  fstream myfile;
  myfile.open(ime_dat,fstream::out);

  std::normal_distribution<double> dist(0,1.0/obj.N);
  std::uniform_int_distribution<int> jpool(0,obj.M-1);
  std::uniform_real_distribution<double> weight(0.0,1.0);
  VectorXd hist = VectorXd::Zero(obj.M);
  VectorXd ksirand = VectorXd::Zero(obj.N);
  for (size_t j = 0; j < obj.nrelax; j++) {
    for (size_t k = 0; k < obj.N; k++) {
      ksirand(k) = dist(gen);
    }
    auto jaj = jpool(gen);
    qmc_step(obj, jaj, weight(gen), ksirand);
    hist +=obj.state;
  }


  for (size_t i = 0; i < obj.nsample; i++) {
    for (size_t k = 0; k < obj.N; k++) {
      ksirand(k) = dist(gen);
    }
    qmc_step(obj, jpool(gen), weight(gen), ksirand);
    if (i%3000 == 0){
      for (size_t i = 0; i < obj.M; i++) {
        myfile << obj.state(i,0) << "\t";
      }
      myfile << "\n";
    }

  }

}



























#endif
