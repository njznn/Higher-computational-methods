#ifndef MD_CMC
#define MD_CMC
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

//POTTS MODEL

template<unsigned int NN>
struct potts {

  int q;
  static const unsigned int N= NN; // size of square lattice

  double J; // exchange integral
  double E=0.0;
  double Esq = 0.0;
  double beta; // beta is 1/T !!
  complex<double> M = 0.0;
  double Msq = 0.0;

  MatrixXd state = MatrixXd::Zero(NN, NN);

  potts(int q_, double J_, double beta_, double seed_): q(q_),J(J_),
   beta(beta_){
     gen.seed(seed_);
     uniform_int_distribution<> dist(1, q);
     for (size_t i = 0; i < N; i++) {
       for (size_t j = 0; j < N; j++) {
         state(i,j) = dist(gen);
       }
     }
     //calculate energy and magnetization of random state
     for (size_t i = 0; i < N; i++) {
       for (size_t j = 0; j < N; j++) {
         E -= J* (((state(i,j)== state((i+1)%N, j))? 1.0:0)
         +((state(i,j)== state(((i-1)==-1? N-1 : i-1), j))? 1.0:0)
         +((state(i,j)== state(i, ((j-1)==-1? N-1 : j-1)))? 1.0:0)
         +((state(i,j)== state(i, (j+1)%N))? 1.0:0));
       }
     }
     E = E*0.5;
     Esq = pow(E, 2);

     for (size_t i = 0; i <N; i++) {
       for (size_t j = 0; j < N; j++) {
          M+= exp((2*numbers::pi * 1i *(state(i,j)-1.0))/((double) q));
       }
     }

     Msq = pow(abs(M), 2);
   }
};

template<unsigned int NN>
double energy_one(potts<NN> & obj, int &i, int &j ){
  return(-obj.J* (((obj.state(i,j)== obj.state((i+1)%obj.N, j))? 1.0:0)
  +((obj.state(i,j)== obj.state(((i-1)==-1? obj.N-1 : i-1), j))? 1.0:0)
  +((obj.state(i,j)== obj.state(i, ((j-1)==-1? obj.N-1 : j-1)))? 1.0:0)
  +((obj.state(i,j)== obj.state(i, (j+1)%obj.N))? 1.0:0)));
}


template<unsigned int NN>
void calc_energy(potts<NN> &obj){
    obj.E = 0.0;
    for (size_t i = 0; i < obj.N; i++) {
      for (size_t j = 0; j < obj.N; j++) {
        obj.E -= obj.J* (((obj.state(i,j)== obj.state((i+1)%obj.N, j))? 1.0:0)
        +((obj.state(i,j)== obj.state(((i-1)==-1? obj.N-1 : i-1), j))? 1.0:0)
        +((obj.state(i,j)== obj.state(i, ((j-1)==-1? obj.N-1 : j-1)))? 1.0:0)
        +((obj.state(i,j)== obj.state(i, (j+1)%obj.N))? 1.0:0));
      }
    }
    obj.E = obj.E*0.5;
}

template<unsigned int NN>
void calc_mag(potts<NN> &obj){
  complex<double> MM = 0;
  for (size_t i = 0; i < obj.N; i++) {
    for (size_t j = 0; j < obj.N; j++) {

      MM += exp((2*numbers::pi * 1i *(obj.state(i,j)-1.0))/((double) obj.q));
      std::cout << exp((2*numbers::pi * 1i *(obj.state(i,j)-1.0))/((double) obj.q)) << '\n';
    }
  }
  obj.M = abs(MM);
}

template<unsigned int NN>
void mc_onemove(potts<NN> &obj, int &ri, int &rj , int &new_q, double &ksi){
  double E_state = energy_one(obj, ri, rj);

  auto old_q = obj.state(ri,rj);
  obj.state(ri, rj) = new_q;

  double E_state_new = energy_one(obj, ri, rj);
  double dE = E_state_new-E_state;

  complex<double> dM = exp((2*numbers::pi * 1i *(new_q-1.0))/((double) obj.q)) -
  exp((2*numbers::pi * 1i *(old_q-1.0))/((double) obj.q));


  if (dE<0){
    obj.E += dE;
    obj.M += dM;
  }
  else if (ksi < exp(-obj.beta*dE)){
    obj.E += dE;
    obj.M += dM;
  }
  else {
    obj.state(ri,rj) = old_q;
  }

}

template<unsigned int NN>
void N_sq_flips(potts<NN> &obj){
  uniform_int_distribution<> dist_ij(0, obj.N-1);
  uniform_real_distribution<> ksi(0.0, 1.0);
  uniform_int_distribution<> qnew(1, obj.q);
  for (size_t i = 0; i < 10000; i++) {
    int q_st = qnew(gen);
    int ri = dist_ij(gen);
    int rj = dist_ij(gen);
    double ks = ksi(gen);
    while(obj.state(ri, rj)==q_st){
      q_st = qnew(gen);
    }
    mc_onemove(obj, ri, rj, q_st, ks);
  }
}

template<unsigned int NN> // N_relax * N^2 steps!
void relax(potts<NN> &obj, int N_relax){
  for (size_t i = 0; i < N_relax; i++) {
    N_sq_flips(obj);
  }
}

template<unsigned int NN>
void potts_simulate_and_flush_beta(int q, double J,double beta_st, int beta_steps, double step_size,
  int Nrelax, int Nsample, std::string ime_dat){
    fstream myfile;
    myfile.open(ime_dat,fstream::out);
    myfile<<"N_relax:"<<Nrelax << "\t Nsample:"<<Nsample<< "\n";
    double beta = beta_st;
    for (size_t i = 0; i < beta_steps; i++) {
      potts<NN> pot(q,J,beta,12);
      relax(pot, Nrelax);
      double avgmag = 0;
      double avgen = 0;
      double avgmagsq = 0;
      double avgensq =0;
      for (size_t j = 0; j < Nsample; j++) {
          N_sq_flips(pot);
          avgmag += abs(pot.M);
          avgen += pot.E;
          avgensq += pow(pot.E, 2);
          avgmagsq += pow(abs(pot.M), 2);
        }

      avgmag = avgmag/(Nsample);
      avgen = avgen/(Nsample);
      avgmagsq = avgmagsq/(Nsample);
      avgensq = avgensq/(Nsample);
      myfile << beta << "\t"<< avgen << "\t" << avgmag << "\t" << avgensq << "\t" << avgmagsq << "\n";
      beta += step_size;

    }
  }

template<unsigned int NN>
void flush_state(potts<NN> &obj, std::string ime_dat){
  fstream myfile;
  myfile.open(ime_dat,fstream::out);
  for (size_t i = 0; i < obj.N; i++) {
    for (size_t j = 0; j < obj.N; j++) {
      myfile << obj.state(i,j) << "\t";
    }
    myfile << "\n";
  }
}





// THIS IS HEISENBERG MODEL WITH CONSTRAINT M=0  for all times!
// N mora biti sodo da lahko upoštevam constraint.
// Initial state will spins in counter clockwise order in x-z plane (z,-x,-z,x)
template<unsigned int NN>
struct heis{
  static const unsigned int N= NN;
  Matrix<double, 3,N> state;
  double Hz =  0.0;
  double J;
  double beta;
  double E=0.0;
  double Esq=0.0;
  Vector3d mag;

  heis(double J_, double Hz_,double beta_ ,double seed_): J(J_), Hz(Hz_), beta(beta_) {
    gen.seed(seed_);
    if (NN>=4 && NN%2==0){
      // initial state
      Vector3d vz(0,0,1);
      Vector3d vx(1,0,0);
      Vector3d vx_(-1,0,0);
      Vector3d vz_(0,0,-1);
      Vector3d arr[4] = {vz, vx_, vz_, vx};
      for (size_t i = 0; i < N; i++) {
        state(all, i) = arr[i%4];
      }
    }
    else {
      std::cout << "Premajhna veriga ali ni sodo število členov" << '\n';
      exit(1);
    }
  }

};

// dont really need it!
template<unsigned int NN>
void calc_mag(heis<NN> &obj){
  obj.mag << 0,0,0;
  for (size_t i = 0; i < obj.N; i++) {
    obj.mag += obj.state(all, i);
  }
  obj.mag = obj.mag/(obj.N);
}

template<unsigned int NN>
void calc_energy(heis<NN> & obj){
  obj.E = 0;
  for (size_t i = 0; i < obj.N; i++) {
    obj.E += -obj.J*(obj.state(all, i).dot(obj.state(all, (i+1)%obj.N))) +
    obj.Hz*obj.state(2, i);
  }
}

Vector3d qrot_vec(Vector3d axis_vec, double theta, const Vector3d &vec){
  if (isnan(theta) || theta ==0){
    return(vec);
  }
  else {
    axis_vec.normalize();
    Quaterniond rot;
    rot.vec() = axis_vec*sin(theta/2);
    rot.w()= cos(theta/2);
    rot.normalize();
    Quaterniond p;
    p.vec() = vec;
    p.w() = 0.0;
    Quaterniond rotatedP = rot* p * (rot.inverse());
    return(rotatedP.vec());
  }
}

template<unsigned int NN>
void mc_step(heis<NN> &obj, int &ri, double &theta, double &ksi){
  Vector3d sum = obj.state(all, ri) + obj.state(all, (ri+1)%obj.N);

  Vector3d Rsi = qrot_vec(sum, theta, obj.state(all, ri));
  Vector3d Rsi_plus = sum - Rsi;

  Vector3d sigi = obj.state(all, ri);
  Vector3d sigi_plus = obj.state(all, (ri+1)%obj.N);

  double dE = -obj.J*(obj.state(all, ((ri-1)==-1? obj.N-1 : ri-1)).dot((Rsi-sigi)) +
  Rsi.dot(Rsi_plus)- sigi.dot(sigi_plus) +
  obj.state(all, (ri+2)%obj.N).dot((Rsi_plus-sigi_plus))) +
  obj.Hz*(Rsi(2) - sigi(2) + Rsi_plus(2)-sigi_plus(2));



  if (dE<0){
    obj.state(all, ri) = Rsi;
    obj.state(all, (ri+1)%obj.N) = Rsi_plus;
    obj.E += dE;
  }
  else if (ksi < exp(-obj.beta*dE)){
    obj.state(all, ri) = Rsi;
    obj.state(all, (ri+1)%obj.N) = Rsi_plus;
    obj.E += dE;
  }
  else {}
}

template<unsigned int NN>
void N_two_flips(heis<NN> &obj){
  uniform_int_distribution<> dist_ij(0, obj.N-1);
  uniform_real_distribution<> ksi(0.0, 1.0);
  uniform_real_distribution<> x(0.0, 1.0);
  for (size_t j = 0; j < obj.N*10; j++) {
      int ri = dist_ij(gen);
      double ks = ksi(gen);
      double theta = 2*numbers::pi*x(gen);
      mc_step(obj, ri, theta, ks);
  }
}

template<unsigned int NN> // N_relax * 2N steps!
void relax(heis<NN> &obj, int N_relax){
  for (size_t i = 0; i < N_relax; i++) {
    N_two_flips(obj);
  }
}

template<unsigned int NN>
VectorXd correlation(heis<NN> &obj, int N_relax, int N_sample){
  VectorXd corr((int)obj.N/2);
  relax(obj, N_relax);
  for (size_t i = 0; i < N_sample; i++) {
    corr += obj.state(2, 0) * obj.state(2,seq(0,(int)obj.N/2));
    N_two_flips(obj);
  }
  corr = corr/N_sample;
  return corr;
}

template<unsigned int NN>
void flush_energy_beta(heis<NN> &obj1, heis<NN> &obj2, heis<NN> &obj3,heis<NN> &obj4,
  heis<NN> &obj5, int N_relax, int N_sample, std::string ime_dat){
  fstream myfile;
  myfile.open(ime_dat,fstream::out);
  myfile << "J=1, h=0 " <<obj1.beta << " "<< obj2.beta<<" "<<obj3.beta << " "<< obj4.beta<<" "<<obj5.beta << "\n ";
  for (size_t i = 0; i < N_sample; i++) {
    myfile<< 2*NN*i <<"\t"<< obj1.E <<"\t" << obj2.E<<"\t"<<obj3.E <<"\t"<< obj4.E <<"\t"<< obj5.E << "\n" ;
    N_two_flips(obj1);
    N_two_flips(obj2);
    N_two_flips(obj3);
    N_two_flips(obj4);
    N_two_flips(obj5);
  }
}






#endif
