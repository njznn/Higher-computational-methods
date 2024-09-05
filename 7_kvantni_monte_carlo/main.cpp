using namespace std;
#include <fstream>
#include <bitset>
using std::ofstream;
#include "qmc.hpp"
#include <iomanip>
#include <chrono>
using namespace std::chrono;

double epsdep(double beta){
  return 0.044* pow(beta,0.505 );
}


int main(int argc, char const *argv[]) {

  // M,B,N,eps,lam,nrelax,nsam,seed,potencial
  qmc qq1(1000,100,1,epsdep(100), 1.0, 100000000,1000000,12,"kvar");

  //relax_and_sample(qq1);
  std::cout << qq1.lam << '\n';
  std::cout << qq1.Eavg << '\n';
  std::cout << qq1.Tavg << '\n';
  std::cout << qq1.Vavg << '\n';
  std::cout << qq1.delez << '\n';

  //histogram:
  //std::string ime_dat = "/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/lhohist_Blam1.txt";
  //relax_and_hist(qq1, ime_dat);


  fstream myfile;
  std::string ime_dat = "/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/2_en_lam_M1000.txt";
  myfile.open(ime_dat,fstream::out);


  VectorXd betasmall = VectorXd::LinSpaced(Sequential,7,0.3,0.9);
  VectorXd beta = VectorXd::LinSpaced(Sequential,50,1,50);
  VectorXd betall(57);
  betall << betasmall, beta;
  VectorXd lam = VectorXd::LinSpaced(Sequential,11,0.0,1.0);

  for (size_t i = 0; i < 11; i++) {
    double epsilon = epsdep(100);
    qmc qq1(1000,100,1,epsilon, lam(i), 10000000,1000000,12,"kvar");
    relax_and_sample(qq1);
    myfile << qq1.Eavg<< '\t' << qq1.Tavg << '\t' << qq1.Vavg<< '\n';


  }



  /*


  //delez dependance of Beta and epsilon
  fstream myfile;
  std::string ime_dat = "/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/1_eps_beta_M100.txt";
  myfile.open(ime_dat,fstream::out);

  VectorXd epsilon = VectorXd::LinSpaced(Sequential,30,0.0,2.0).transpose();



  myfile<< 0.1 <<"\t"<<1.0 <<"\t"<<2.0 <<"\t"<<3.0 <<"\t"<<5.0 <<"\t"<<
  10.0 <<"\t"<<20.0 <<"\t"<<30.0 <<"\n";
  //myfile<< 50 <<"\t"<<100 <<"\t"<<200 <<"\t"<<400 <<"\t"<<600 <<"\t"<<
  //800 <<"\t"<<1000 <<"\t"<<1200 <<"\n";
  //myfile<< 0.01 <<"\t"<<0.1 <<"\t"<<0.2 <<"\t"<<0.4 <<"\t"<<0.6 <<"\t"<<
  //0.8 <<"\t"<<1.0 <<"\t"<<2.0 <<"\n";

  for (size_t i = 0; i < 30; i++) {
    std::cout << epsilon(i) << '\n';

    qmc qq1(100,0.1,1,epsilon(i), 0, 1000000,1000000,12,"kvar");
    qmc qq2(100,1.0,1,epsilon(i),0, 1000000,1000000,12,"kvar");
    qmc qq3(100,2.0,1,epsilon(i), 0, 1000000,1000000,12,"kvar");
    qmc qq4(100,3.0,1,epsilon(i), 0, 1000000,1000000,12,"kvar");
    qmc qq5(100,5.0,1,epsilon(i), 0, 1000000,1000000,12,"kvar");
    qmc qq6(100,10.0,1,epsilon(i), 0, 1000000,1000000,12,"kvar");
    qmc qq7(100,20.0,1,epsilon(i), 0, 1000000,1000000,12,"kvar");
    qmc qq8(100,30.0,1,epsilon(i), 0, 1000000,1000000,12,"kvar");

    relax_and_sample(qq1);
    relax_and_sample(qq2);
    relax_and_sample(qq3);
    relax_and_sample(qq4);
    relax_and_sample(qq5);
    relax_and_sample(qq6);
    relax_and_sample(qq7);
    relax_and_sample(qq8);

    myfile<<epsilon(i) <<"\t"<<qq1.delez <<"\t"<<qq2.delez <<"\t"<<qq3.delez <<"\t"<<qq4.delez <<"\t"<<qq5.delez <<"\t"<<
    qq6.delez <<"\t"<<qq7.delez <<"\t"<<qq8.delez <<"\n";

  }

  */




  return 0;
}
