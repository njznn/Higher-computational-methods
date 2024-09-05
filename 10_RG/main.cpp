using namespace std;
#include <fstream>
using std::ofstream;
#include "rg.hpp"
#include <iomanip>
#include <chrono>
#include <stdlib.h>
using namespace honeycomb;
/*
void output_numpy_arrays_to_txt_file(const VectorXd& array1,
  const VectorXd& array2, std::string filename) {
  std::ofstream f(filename);
  for (int i = 0; i < array1.size(); ++i) {
    f << array1(i) << " " << array2(i) << "\n";
  }
  f.close();
}


void print_eigen_matrix_to_txt_file(const Eigen::MatrixXd& matrix, const string& filename) {
  ofstream outfile(filename);
  for (int i = 0; i < matrix.rows(); ++i) {
    for (int j = 0; j < matrix.cols(); ++j) {
      outfile << matrix(i, j) << " ";
    }
    outfile << endl;
  }
  outfile.close();
}


void Niterdep() {
  VectorXd koraki(4);
  koraki<<2,4,10,20;
  VectorXd betasm = VectorXd::LinSpaced(39, 0.05, 1.0);
  VectorXd betalg = VectorXd::LinSpaced(120, 1.025, 4.0);
  VectorXd beta(betasm.size()+ betalg.size());
  beta <<betasm , betalg;

  MatrixXd allinf(159,5);
  allinf(all, 0) = beta;

  for (int i = 0; i < koraki.size(); ++i) {
    allinf(all, i+1) = lnZ_sweep_beta(6,1,koraki(i),20);


  }
  print_eigen_matrix_to_txt_file(allinf, "Nitdep_q6_J05_M20.txt" );
}

void Mdep(){
  VectorXi M(6);
  M<<2,4,6,10,16,20;
  VectorXd betasm = VectorXd::LinSpaced(39, 0.05, 1.0);
  VectorXd betalg = VectorXd::LinSpaced(120, 1.025, 4.0);
  VectorXd beta(betasm.size()+ betalg.size());
  beta <<betasm , betalg;

  MatrixXd allinf(159,7);
  allinf(all, 0) = beta;

  for (int i = 0; i < M.size(); i++) {
    allinf(all, i+1) = lnZ_sweep_beta(6,1,3,M(i));
    std::cout << M(i) << '\n';


  }
  print_eigen_matrix_to_txt_file(allinf, "Mdep_q6_J05_Nit3.txt" );
}

*/

int main() {
  //Partitsum(2,1,1,1,2);

  lnZ_sweep_beta_tofile(5,2,6,10,"heks__q=5_N12_M26_J1_nonorm.txt");
  //Niterdep();
  // Print the duration in milliseconds.
  VectorXd betasm = VectorXd::LinSpaced(39, 0.05, 1.0);
  VectorXd betalg = VectorXd::LinSpaced(120, 1.025, 4.0);
  VectorXd beta(betasm.size()+ betalg.size());
  beta <<betasm , betalg;
  //std::cout << beta << '\n';
  //VectorXd res = lnZ_sweep_beta(4,1,16,18);
  //output_numpy_arrays_to_txt_file(beta, res, "lnz_12st_q=4_M20_cpp.txt");
  //Mdep();
  //std::cout << Partitsum(6,1,1,16,25) << '\n';
  return 0;
}
