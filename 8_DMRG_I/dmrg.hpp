#ifndef DMRG_HPP
#define DMRG_HPP
#include <unsupported/Eigen/KroneckerProduct>
#include <spectra/include/Spectra/SymEigsSolver.h>
#include <spectra/include/Spectra/MatOp/SparseGenMatProd.h>
using namespace Spectra;
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include <iostream>
#include <bitset>
using namespace Spectra;
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <spectra/include/Spectra/GenEigsSolver.h>





/// CONSTRUCTING HAMILTONIAN:

//one way:
/*Matrix4d const h;
h << 1,0,0,0,
     0,-1,2,0,
     0,2,-1,0,
     0,0,0,1;
*/
MatrixXd hij( Matrix4d &h, int j, int  n){
    return kroneckerProduct(
      kroneckerProduct(MatrixXd::Identity(pow(2, (j-1)),pow(2, (j-1))),h).eval(),
      MatrixXd::Identity(pow(2, (n-j-1)), pow(2, (n-j-1)))).eval();
 }

MatrixXd Heisenberg_H(int n, std::string robni){
  Matrix4d h;
  h << 1,0,0,0,
       0,-1,2,0,
       0,2,-1,0,
       0,0,0,1;

  MatrixXd Hmat((int) pow(2, n),(int)  pow(2,n));
  for (size_t i = 1; i < n; i++) {
    Hmat += hij(h, i, n);
  }
  if (robni=="periodic"){
    Matrix2d Sx;
    Matrix2d Sz;
    Matrix2cd Sy;
    Sx << 0,1,1,0;
    Sz << 1,0,0,-1;
    Sy << 0,-1i,1i,0;

    Hmat += kroneckerProduct(Sx,kroneckerProduct(MatrixXd::Identity(pow(2, (n-2)),pow(2, (n-2))), Sx).eval()).eval()
    +kroneckerProduct(Sz,kroneckerProduct(MatrixXd::Identity(pow(2, (n-2)),pow(2, (n-2))), Sz).eval()).eval()
    +kroneckerProduct(Sy,kroneckerProduct(MatrixXd::Identity(pow(2, (n-2)),pow(2, (n-2))), Sy).eval()).eval().real();

  }
  return Hmat;
}

VectorXd basestate(MatrixXd & M){
  DenseSymMatProd<double> op(M);
  SymEigsSolver<DenseSymMatProd<double>> eigs(op, 1, 6);
  eigs.init();
  int nconv = eigs.compute(SortRule::SmallestAlge );
  if (eigs.info() == CompInfo::Successful){

    return eigs.eigenvectors();
  }
  else{
    std::cout << "error" << '\n';
    exit(1);
  }
}


// drugi nacin
template <int T>
Vector<complex<double>, Dynamic> H_on_state(VectorXcd state , std::string robni){
  static constexpr int N = T;
  Vector<complex<double>, Dynamic> newvec = Matrix<complex<double>, Dynamic, 1> ::Zero(pow(2, N),1);
  int dolz;
  long int nstate = pow(2, N);
  dolz = (robni != "periodic") ? N-1 : N;
  for (size_t i = 0; i < dolz; i++) {
    for (unsigned int ind= 0; ind < nstate; ind++) {
      auto bind = bitset<N>(ind);
      string twobits = (to_string(bind[(i+1)% N]) + to_string(bind[i]));
      if (twobits == "11" || twobits=="00"){
        newvec(ind) +=  state(ind);
      }
      else if (twobits == "01"){
        auto swapped = bind.flip(i).flip((i+1)%N);
        int swapind = stoi(swapped.to_string(), nullptr, 2);
        newvec(ind) += (-1.0)*state(ind)
        + 2.0* state(swapind);

      }

      else if(twobits == "10"){
        auto swapped = bind.flip(i).flip((i+1)%N);
        int swapind = stoi(swapped.to_string(), nullptr, 2);
        newvec(ind) += (-1.0)*state(ind)
        + 2.0* state(swapind);
      }
    }

  }
  return newvec;
}


VectorXi getnonzeroind(VectorXd  vec){
  VectorXi ind(vec.size());
  int mesto = 0;

  for (size_t i = 0; i < vec.size(); i++) {
    if (abs(vec(i)) > pow(10, -6)){
      ind(mesto) = i;
      mesto +=1;
    }
  }

  return ind(seq(0,mesto-1));
}


template <int T>
MatrixXd makeHeisenH(std::string robni){
  int N = pow(2, T);
  MatrixXd Hmat(N,N);
  VectorXd estate = Matrix<double, Dynamic, 1> ::Zero(N,1);
  for (size_t i = 0; i < N; i++) {
    estate(i)=1;
    if (i>0){
      estate(i-1) = 0;
    }
    Hmat(i, all) = H_on_state<T>(estate, robni).real();
  }

  return Hmat;

}

template <int T>
SparseMatrix<double> makeSparseH(std::string robni){
  const int N = pow(2, T);
  SparseMatrix<double> M(N,N);
  std::vector<Triplet<double>> tripletList;
  VectorXd estate = Matrix<double, Dynamic, 1> ::Zero(N,1);
  for (size_t i = 0; i < N; i++) {
    estate(i)=1;
    if (i>0){
      estate(i-1) = 0;
    }
    VectorXd vec = H_on_state<T>(estate, robni).real();
    VectorXi ind = getnonzeroind(vec);

    for (size_t j = 0; j < ind.size(); j++) {
      tripletList.push_back(Triplet<double>(i,ind(j),vec(ind(j))));
    }
  }
  M.setFromTriplets(tripletList.begin(), tripletList.end());

  return M;
}

VectorXd basestatefromsparse(SparseMatrix<double> & M){
  SparseGenMatProd<double> op(M);
  SymEigsSolver<SparseGenMatProd<double>> eigs(op, 1, 6);
  eigs.init();
  int nconv = eigs.compute(SortRule::SmallestAlge);
  if (eigs.info() == CompInfo::Successful){
    return eigs.eigenvectors();
  }
  else{
    std::cout << "error" << '\n';
    exit(1);
  }
}








#endif
