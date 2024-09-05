using namespace std;
#include <fstream>
#include <bitset>
using std::ofstream;
#include "cmc.hpp"
#include <iomanip>
#include <chrono>
using namespace std::chrono;
#include <omp.h>
#include <Eigen/Geometry>




int main() {
  //std::string ime_dat = "/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_32_q=2_J=05.txt";

  //potts_simulate_and_flush_beta<32>(3,0.5,0.01,40,0.1,
  //10000, 500, ime_dat);



  //fstream myfile;
  //myfile.open(ime_dat,fstream::out);
  //for (size_t i = 0; i < 1000; i++) {
    /* code */
  //}

  //N_sq_flips(pot);




  // Heisenberg:

  std::string ime_dat = "/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/heis_state100_b1000.txt";
  fstream myfile;
  myfile.open(ime_dat,fstream::out);
  //relax dependance:
  //J, Hz, beta, seed

  /*
  heis<100> hai2(1,0,0.1,12);
  heis<100> hai3(1,0,1,12);
  heis<100> hai4(1,0,10,12);
  heis<100> hai5(1,0,20,12);
  */
  //std::cout << correlation(hai, 10000, 10000) << '\n';
  //flush_energy_beta(hai1, hai2, hai3, hai4, hai5, 0, 1000000, ime_dat);
  double h = 0.0;
  heis<100> hai1(1,h,1000,12);
  auto v1 = correlation(hai1, 100000, 10000);
  /*
  heis<500> hai2(1,h,0.1,12);
  auto v2 = correlation(hai2, 100000, 10000);
  heis<500> hai3(1,h,1,12);
  auto v3 = correlation(hai3, 100000, 10000);
  heis<500> hai4(1,h,2,12);
  auto v4 = correlation(hai4, 100000, 10000);
  heis<500> hai5(1,h,6,12);
  auto v5 = correlation(hai5, 100000, 10000);
  heis<500> hai6(1,h,30,12);
  auto v6 = correlation(hai6, 100000, 10000);

  */
  for (size_t i = 0; i < 100; i++) {
      myfile<< hai1.state(0, i)<<"\t"<<hai1.state(1,i)<<"\t"
      << hai1.state(2, i)<< "\n";
  }

/*
  myfile << "J=1, beta=0, 10^5" << "10^6"<< "..." << "10^8 \n";
  for (size_t i = 0; i < hai1.N/2; i++) {
    myfile<< v1(i) <<"\t"<< v2(i) <<"\t" << v3(i)<<"\t"<<v4(i)
    <<"\t"<< v5(i)<<"\t"<<v6(i)<< "\n" ;

  }
*/


  return 0;
}
