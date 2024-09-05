#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include <math.h>
#include <complex.h>
#include <fstream>
#include <bitset>
using std::ofstream;
#include "qb.hpp"
#include <iomanip>



int main() {

  qb arr[] = {qb(0.01, "cas"), qb(0.01, "cas"), qb(0.01, "cas"),
              qb(0.01, "cas"), qb(0.01, "cas"), qb(0.01, "cas")};


  //qb chain(0.1, "cas");


  //propagate_time(chain, "S1",200);

  //std::cout <<  Zb<6>(arr, 5, 1)<< '\n';
  //write_mat_txt(arr[0],Zb<6>(arr, 5, 0),
  // "/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/Z(B)_N=12_all_real.txt" );
  //std::cout <<  sigmaz_autocol_one(chain, 0, 2) << '\n';
  //write_mat_txt(arr[0],Cz<6>(arr, 2, 5),
  //"/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C22(t)_N20_S2_01_r.txt" );

  write_mat_txt(arr[0],spin_current<6>(arr, 50),
  "/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/J(t)_N12_S2_01_r_long.txt" );

  return 0;
}
