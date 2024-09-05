using namespace std;
#include <fstream>
using std::ofstream;
#define SIMPLE_LINEAR_REGRESSION_IMPLEMENTATION
#include "simple_linear_regression.h"
#include "tebd.hpp"
#include <iomanip>
#include <chrono>
#include <stdlib.h>

typedef std::vector<Matrix<complex<double>, Dynamic, Dynamic>> vecmat;

/*
double calc_afm_ground_energy(vecmat &Bijinit,
   vecmat &laminit, VectorXi randmstate_neel, int Nsteps, double stepsize, double &err){
     auto propbij = Bijinit;
     auto proplam = laminit;

     int Nspin = randmstate_neel.size();

     suz_trotter(propbij, proplam, pow(10, -7), err,(-1.0)*stepsize,
     Nspin,Nsteps, "S2");

     double y[10];
     double x[10];
     double betafin = Nsteps*stepsize;
     double slope=0.0;
     for (size_t i = 1; i < 11; i++) {
       suz_trotter(propbij, proplam, pow(10, -8), err,-0.01 ,Nspin,10, "S2");
       y[i-1] = log(abs(coefofstateB<complex<double>>(propbij, randmstate_neel, proplam, Nspin, 2))/(
       abs(coefofstateB<complex<double>>(Bijinit, randmstate_neel, laminit, Nspin,2))));
       x[i-1] = betafin + 0.01*i;
     }


     simple_linear_regression(x, y, sizeof(x) / sizeof(x[0]), &slope, NULL, NULL, NULL, NULL, NULL);

     return slope/(10*Nspin);
   }



double calc_afm_encur(vecmat &Bijinit,
   vecmat &laminit, VectorXi randmstate_neel, double &err, complex<double> beta){
     auto propbij = Bijinit;
     auto proplam = laminit;

     int Nspin = randmstate_neel.size();
     complex<double> notrel = 0.0;
     //suz_trotter(propbij, proplam, pow(10, -7), err,(-1.0)*stepsize,
     //Nspin,Nsteps, "S2", curz);

     double y[10];
     double x[10];
     double slope=0.0;
     for (size_t i = 1; i < 11; i++) {
       suz_trotter(propbij, proplam, pow(10, -8), err,-0.01 ,Nspin,10, "S2", notrel);
       y[i-1] = log(abs(coefofstateB<complex<double>>(propbij, randmstate_neel, proplam, Nspin, 2))/(
       abs(coefofstateB<complex<double>>(Bijinit, randmstate_neel, laminit, Nspin,2))));
       x[i-1] = beta.real() + 0.01*i;
     }


     simple_linear_regression(x, y, sizeof(x) / sizeof(x[0]), &slope, NULL, NULL, NULL, NULL, NULL);

     return slope/(10);
   }

/*
VectorXd calc_n_afm_en_dep(){
  VectorXi n(1);
  VectorXd energije(1);
  n<< 100;
  double err = 0;
  for (size_t j = 0; j <n.size(); j++) {
    VectorXi neel(n(j));
    for (size_t i = 0; i < neel.size(); i++) {
      neel(i) = i%2;
    }
    std::cout << n(j) << '\n';
    auto res = MPSrandomstate<complex<double>>(n(j), 128, "False");
    double en = calc_afm_ground_energy(res[0], res[1],neel, 10,0.1,err);
    energije(j) = en;
  }
  return energije;
}


MatrixXd betadepn(){
  VectorXi n(10);
  n << 4,6,8,10,12,14,16,18,20,32;
  VectorXd beta(26);
  VectorXd betabig = VectorXd::LinSpaced(17,1,9);
  VectorXd betasm = VectorXd::LinSpaced(9,0.1,0.9);
  beta << betasm, betabig;
  MatrixXd dep(int(beta.size()), int(n.size() + 1));
  dep(all, 0) = beta;
  for (size_t i = 0; i < n.size(); i++) {
    double err = 0;
    VectorXi neel(n(i));
    for (size_t j = 0; j < neel.size(); j++) {
      neel(j) = j%2;
    }
    auto res = MPSrandomstate<complex<double>>(n(i), 128, "False");
    complex<double> curz = 0.0;
    int ind = 0;
    for (size_t k = 0; k < 90; k++) {
      suz_trotter<complex<double>>(res[0], res[1], pow(10, -8), err,-0.1 ,n(i),1, "S2", curz);
      if ((curz.real()- beta(ind))< pow(10, -5)){
        dep(ind, i+1)= calc_afm_encur(res[0], res[1], neel, err, curz);
        ind +=1;
      }
    }

  }

  return dep;
}


MatrixXd M_e0_dep(){
  VectorXi Nspin(3);
  Nspin <<14;
  MatrixXd en(6,3);
  VectorXd size(6);
  for (size_t i = 0; i < Nspin.size(); i++) {
    VectorXi neel(Nspin(i));
    for (size_t j = 0; j < neel.size(); j++) {
      neel(j) = j%2;
    }
    double err = 0;
    size << 16,32,64,128,256,512;
    std::cout << i << '\n';
    for (size_t j = 0; j < size.size(); j++) {
      auto res = MPSrandomstate<complex<double>>(Nspin(i), size(j), "False");
      en(j,i) = calc_afm_ground_energy(res[0],res[1], neel, 50, 0.1, err);
    }
  }
  return en;
}

MatrixXd calc_corr(int Nsp){
  double err = 0.0;
  int Nspin = Nsp;
  std::cout << "beta1" << '\n';
  std::cout << "N" << Nsp << '\n';
  std::cout << "M=80 " << '\n';
  MatrixXd corr(Nspin, 1);
  VectorXi mesta(1);
  mesta << Nspin/2;
  for (size_t j = 0; j < mesta.size(); j++) {
    std::cout << j << '\n';
    auto res = MPSrandomstate<complex<double>>(Nspin, 50, "False");
    suz_trotter(res[0], res[1], pow(10, -7), err,-0.1 ,Nspin,10, "S2");
    for (size_t i = 0; i < Nspin; i++) {
      std::cout << i << '\n';
      corr(i,j) = sigmaksigmaj<complex<double>>(res[0], res[0], res[1], res[1],2,mesta(j),i);

    }
    corr(all, j) /= corr(all, j).cwiseAbs().maxCoeff();

  }
  return corr;
}

*/
MatrixXd calc_corr__time(int Nsp){
  double err = 0.0;
  int Nspin = Nsp;
  VectorXd time = VectorXd::LinSpaced(30,0,29);
  MatrixXd corr(time.size(), 3);
  VectorXi mesta(3);
  mesta << 0, Nspin/2, Nspin-1;
  std::cout << "Nspin=" <<Nspin/2<< '\n';
  auto res = MPSrandomstate<complex<double>>(Nspin, 50, "False");

  //suz_trotter(res[0], res[1], pow(10, -7), err,-0.1 ,Nspin,20, "S2");
  auto copy = res;
  sigma_zj_MPA(copy[0],0);
  for (size_t k = 0; k < mesta.size(); k++) {
    for (size_t j = 0; j < time.size(); j++) {
      std::cout << j << '\n';
      corr(j,k) = sigmaksigmaj<complex<double>>(res[0], copy[0], res[1], copy[1],2,mesta(k),0);
      suz_trotter(res[0], res[1], pow(10, -7), err,-0.1i ,Nspin,10, "S2");
      suz_trotter(copy[0], copy[1], pow(10, -7), err,-0.1i ,Nspin,10, "S2");
    }
    corr(all, k)/= corr(all, k).cwiseAbs().maxCoeff();
  }

  return corr;
}

/*
MatrixXd domainevol(VectorXd casi){
  VectorXi domainconf(4);
  int Nspin = 60;
  domainconf << 20,10,5,25;
  double err = 0;
  auto res = MPSdomainstate<complex<double>>(Nspin, 50, "False", domainconf, pow(10, -2));
  MatrixXd szprofil(Nspin, casi.size());
  double tcur =0.0;
  int ind = 0;
  while (tcur <= casi(casi.size()-1)) {
    std::cout << tcur << '\n';
    if (abs(tcur - casi(ind)) < pow(10,-4)){
      for (size_t i = 0; i < Nspin; i++) {
        szprofil(i,ind) = sigmaksigmaj(res[0],res[0], res[1], res[1], 2, i, -1);
      }
      szprofil(all,ind) /= szprofil.maxCoeff();
      ind +=1;
      std::cout << ind << '\n';
    }
    suz_trotter<complex<double>>(res[0], res[1],pow(10, -20),err, 0.1i, Nspin, 1, "S2");
    tcur += 0.1;

  }

  return szprofil;
}

/*
  MatrixXd normerr(){
    VectorXi M(5);
    int Nspin = 20;
    std::cout << "Nspin32" << '\n';
    VectorXd time = VectorXd::LinSpaced(20,0,19);
    MatrixXd sol(int(time.size()), int(M.size()));
    M<< 32,50,64,80,100;
    for (size_t i = 1; i <2; i++) {
      std::cout << "M:" <<M(i)<< '\n';
      VectorXi domainconf(2);
      domainconf << 5,5;
      double err = 0;
      auto res = MPSrandomstate<complex<double>>(Nspin, M(i), "False");
      double initialnorm = TscalprodB<complex<double>>(res[0], res[0], res[1], res[1], 2);
      std::cout << initialnorm << '\n';
      sol(0,i) = 0;
      for (size_t j = 1; j < 20; j++) {
        suz_trotter<complex<double>>(res[0], res[1],pow(10, -30),err, 0.1i, Nspin, 10, "S2");
        double norm = TscalprodB<complex<double>>(res[0], res[0], res[1], res[1], 2);
        sol(j,i) = norm/initialnorm-1;
        std::cout << initialnorm << '\n';
        std::cout << "----" << '\n';
        std::cout << norm << '\n';
      }

    }
    return sol;
  }

  template<int Nsp>
  VectorXd ksidep(){
    VectorXd ksi(8);
    ksi<<0.0, pow(10,-7), pow(10,-6), pow(10,-5), pow(10,-4), pow(10,-3), pow(10,-2),
    pow(10,-1);
    VectorXd sol(ksi.size());


    VectorXi domainconf(2);
    domainconf << Nsp/2,Nsp/2;
    srand(12);
    auto res = MPSdomainstate<complex<double>>(Nsp, 50, "False", domainconf, 0.0);
    auto koefzero = veckoeffrmpsB<complex<double>, Nsp>(res[0], res[1], 2);
    for (size_t k = 0; k < ksi.size(); k++) {
      srand(12);
      auto res2 = MPSdomainstate<complex<double>>(Nsp, 50, "False", domainconf, ksi(k));
      auto koef = veckoeffrmpsB<complex<double>, Nsp>(res2[0], res2[1], 2);
      sol(k) = pow((koef-koefzero).cwiseAbs().maxCoeff(),2) / pow(koefzero.cwiseAbs().maxCoeff(), 2);
    }

    return sol;
  }
  */
int main(int argc, char const *argv[]) {
  /*
  int N10 = pow(2,4);
  VectorXcd base10(N10);
  base10.real() = VectorXd::Random(N10);
  //base10.imag() = VectorXd::Random(N10);
  base10.normalize();
  VectorXcd base10sec(N10);
  base10sec.real() = VectorXd::Random(N10);
  base10sec.imag() = VectorXd::Random(N10);
  base10sec.normalize();
  VectorXcd conj = base10.conjugate();

  std::cout << base10 << '\n';
  std::cout << "---" << '\n';
  //std::cout << calc_corr(100) << '\n';
  auto mpsB = MPSwB<complex<double>>(base10, 4, 2, "ne", 0);
  sigma_zj_MPA<complex<double>>(mpsB[0], 2);
  std::cout << veckoeffrmpsB<complex<double>,4>(mpsB[0], mpsB[1], 2) << '\n';



  //suz_trotter(res[0], res[1], pow(10, -7), err,-0.1 ,10,100, "S2");
  */
  std::cout << "korelacija 32 se enkrat" << '\n';
  std::cout << calc_corr__time(32) << '\n';
  std::cout << "----" << '\n';
  /*
  VectorXd casi(18);
  casi << 0.0,0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10,11,12,13,14,15,16;//17,18,19,20;//,21,22,23,24,25;
  std::cout << casi << '\n';
  std::cout << "3 domains,60, diff size!" << '\n';
  std::cout << domainevol(casi) << '\n';
  //std::cout << normerr() << '\n';
  /*
  MatrixXd sol(8,7);
  std::cout << "start" << '\n';
  sol(all, 0) = ksidep<8>();
  sol(all, 1) = ksidep<10>();
  sol(all, 2) = ksidep<12>();
  sol(all, 3) = ksidep<14>();
  sol(all, 4) = ksidep<16>();
  std::cout << "18" << '\n';
  sol(all, 5) = ksidep<18>();
  sol(all, 6) = ksidep<20>();

  std::cout << sol << '\n';
  std::cout << "break" << '\n';


  //auto res2 = MPSrandomstate<complex<double>>(4, 64, "False");


  //std::cout << veckoeffrmpsB<complex<double>,6>(res[0], res[1],2 ) << '\n';
  /*
  double err = 0;
  VectorXi domainconf(2);
  domainconf << 20,20;


  auto res = MPSdomainstate<complex<double>>(40, 64, "False", domainconf);
  VectorXd szprofil(40);
  suz_trotter<complex<double>>(res[0], res[1],pow(10, -20),err, 0.1i, 40, 50, "S2");
  std::cout << "----" << '\n';
  for (size_t i = 0; i < 40; i++) {
    szprofil(i) = sigmaksigmaj(res[0],res[0], res[1], res[1], 2, i, -1);
  }
  szprofil /= szprofil.maxCoeff();
  std::cout << szprofil << '\n';


  //suz_trotter<complex<double>>(res[0], res[1], pow(10, -7), err,0.1i ,4,1, "S2");
  std::cout << "break" << '\n';
  //std::cout << calc_corr__time(32) << '\n';



  //std::cout << res[0][9] << '\n';

  //std::cout << res[0][7] << '\n';
  //std::cout << res[1][7] << '\n';

  //std::cout << betadepn() << '\n';
  //auto vec1 = veckoeffrmps<complex<double>, 10>(resA[0], 2);
  //std::cout << vec1 << '\n';

  //auto vec2 = veckoeffrmpsB<complex<double>, 10>(res[0],res[1],2);

  //std::cout << vec2 << '\n';



  //auto vec1 = veckoeffrmpsB<complex<double>, 10>(res[0],res[1], 2);
  //std::cout << TscalprodB<complex<double>>(copy[0],res[0],copy[1],res[1], 2) << '\n';
  //std::cout << M_e0_dep() << '\n';


  /*
  fstream myfile;
  std::string ime_dat = "/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/9_DMRG_II/randvec.txt";
  myfile.open(ime_dat,fstream::out);


  for (size_t i = 0; i < base10.size(); i++) {
    myfile << base10.real()(i) << '\t' << base10.imag()(i)<< '\t'<< base10sec.real()(i)<< '\t'<<base10sec.imag()(i) << '\n';
  }
  */

  return 0;
}
