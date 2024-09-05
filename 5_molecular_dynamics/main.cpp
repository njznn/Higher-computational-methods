using namespace std;
#include <fstream>
#include <bitset>
using std::ofstream;
#include "md.hpp"
#include <iomanip>
#include <random>
#include <chrono>
using namespace std::chrono;
#include <omp.h>

int main() {

  random_device rd{};
  mt19937 gen{rd()};
  auto start = high_resolution_clock::now();
  //N ,dt, TL,TR,lambda,tdead,tmax, tau
  // OPTIMAL PARAMETERS FOR NO0SE-HOVER:
  // dt = 0.1, trelax=tsample=100000, tau = 0.5

  int Num = 10;
  mch ch0(Num,0.1,1.0,2.0,0,100000,100000,0.5);
  mch ch1(Num,0.1,1.0,2.0,0.01,100000,100000,0.5);
  mch ch2(Num,0.1,1.0,2.0,0.1,100000,100000,0.5);
  mch ch3(Num,0.1,1.0,2.0,0.2,100000,100000,0.5);
  mch ch4(Num,0.1,1.0,2.0,0.4,100000,100000,0.5);
  mch ch5(Num,0.1,1.0,2.0,0.6,100000,100000,0.5);
  mch ch6(Num,0.1,1.0,2.0,0.8,100000,100000,0.5);
  mch ch7(Num,0.1,1.0,2.0,1.0,100000,100000,0.5);


  RK4_iteration(ch0);
  RK4_iteration(ch1);
  RK4_iteration(ch2);
  RK4_iteration(ch3);
  RK4_iteration(ch4);
  RK4_iteration(ch5);
  RK4_iteration(ch6);
  RK4_iteration(ch7);



/*
  //OPTIMAL PARAMETERS FOR MAXWELL:
  //dt = .0.1, trelax=tsample=100000, \tau = 0.5

  int Num = 10;
  mch ch0(Num,0.1,1.0,2.0,0,100000,100000,0.5);
  mch ch1(Num,0.1,1.0,2.0,0.01,100000,100000,0.5);
  mch ch2(Num,0.1,1.0,2.0,0.1,100000,100000,0.5);
  mch ch3(Num,0.1,1.0,2.0,0.2,100000,100000,0.5);
  mch ch4(Num,0.1,1.0,2.0,0.4,100000,100000,0.5);
  mch ch5(Num,0.1,1.0,2.0,0.6,100000,100000,0.5);
  mch ch6(Num,0.1,1.0,2.0,0.8,100000,100000,0.5);
  mch ch7(Num,0.1,1.0,2.0,1.0,100000,100000,0.5);


  RK4_maxwell_prop(ch0);
  RK4_maxwell_prop(ch1);
  RK4_maxwell_prop(ch2);
  RK4_maxwell_prop(ch3);
  RK4_maxwell_prop(ch4);
  RK4_maxwell_prop(ch5);
  RK4_maxwell_prop(ch6);
  RK4_maxwell_prop(ch7);


  */
  /*
  string ime_dat = "/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/maxwell_lamda_N70.txt";
  fstream myfile;
  myfile.open(ime_dat,fstream::out);
  myfile << "N = " << ch0.N<<" " <<"dt = "<<ch0.dt <<"lambda ="<<ch0.lambda<<" " <<
  "tdead"<< ch0.tdead <<" " <<"tmax="<< ch0.tmax << " "<< "tau="<< ch0.tau<<  "\n";

  myfile << "lambda = 0.0, 0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0 " << "\n";
  //myfile << "tau = 0.1, 0.5, 1.0, 1.5, 2.0" << "\n";
  for (size_t i = 0; i < ch0.N; i++) {
    myfile << ch0.Tavg(i) << "\t" <<ch0.Javg(i)<< "\t"<< ch1.Tavg(i) <<"\t"<< ch1.Javg(i)
     << "\t" << ch2.Tavg(i) << "\t" <<ch2.Javg(i)<< "\t"<< ch3.Tavg(i) <<"\t" <<
    ch3.Javg(i)<<"\t"<< ch4.Tavg(i) << "\t" <<ch4.Javg(i)<< "\t"<< ch5.Tavg(i) <<"\t"<< ch5.Javg(i)
     << "\t" << ch6.Tavg(i) << "\t" <<ch6.Javg(i)<< "\t"<< ch7.Tavg(i) << "\t"<<
    ch7.Javg(i)<< "\n";
    }
    */
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    std::cout << duration.count() << '\n';
    std::cout << ch7.state << '\n';

  return 0;
}
