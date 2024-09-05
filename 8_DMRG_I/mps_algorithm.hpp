#ifndef MPS_ALGORITHM_HPP
#define MPS_ALGORITHM_HPP
#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include <Eigen/SVD>
#include <cmath>
#include <bitset>
#include <string>
#include <unsupported/Eigen/KroneckerProduct>
#include <complex>



// MPS decomposition
// ACTUAL TRESHOLD VALUE is given as treshold * |max eigenvalue|
// lst value of eigenvalue is zero because no svd is made at last spin,
//only svd V is assigned to An!, same for psiarr
template <typename T>
std::vector<std::vector<Matrix<T, Dynamic, Dynamic>>>
MPSdecomp(Vector<T, Dynamic>  &stateref,int Nspin,
  int d, std::string odrezi, double treshold){

  Matrix<T, Dynamic, Dynamic, RowMajor> state = stateref;

  int mi = 1;

  std::vector<Matrix<T, Dynamic, Dynamic>> matrices(Nspin);
  std::vector<Matrix<T, Dynamic, Dynamic>> psiarr(Nspin-1);
  std::vector<Matrix<T, Dynamic, Dynamic>> singularval(Nspin-1);
  std::vector<Matrix<T, Dynamic, Dynamic>> singularvalcut(Nspin-1);

  int Nstates = (int) pow(d, Nspin);

  for (size_t i = 0; i < Nspin-1; i++) {

    state.resize(d*mi, pow(d, Nspin-1-i) );

    BDCSVD<Matrix<T, Dynamic, Dynamic>> bdcsvd( state,  ComputeThinU | ComputeThinV);

    mi = bdcsvd.singularValues().size();

    if (odrezi=="odrezi"){
      double maxsing = bdcsvd.singularValues().maxCoeff();
      bdcsvd.setThreshold(treshold/abs(maxsing));
      mi = bdcsvd.rank();

    }

    matrices[i] = bdcsvd.matrixU()(all, seq(0,mi-1));

    state.resize(mi, pow(d, Nspin-1-i) );
    for (size_t j = 0; j < mi; j++) {

      state(j, all) = bdcsvd.singularValues()(j)* (bdcsvd.matrixV().adjoint())(j, all);
    }

    psiarr[i] = state;
    singularval[i] = bdcsvd.singularValues();
    singularvalcut[i] = bdcsvd.singularValues()(seq(0,mi-1));
  }

  matrices[Nspin-1] = state.transpose();


  std::vector<std::vector<Matrix<T, Dynamic, Dynamic>>> res = {matrices, psiarr,
    singularval,  singularvalcut};
  return res;


}


template <typename T>
VectorXd entangeled_entropy( std::vector<Matrix<T, Dynamic, Dynamic>>  & singularvalcut){

  VectorXd entropy(size(singularvalcut));

  double Si = 0;
  double truncerr = 0;
  for (size_t i = 0; i < size(singularvalcut); i++) {
    for (size_t j = 0; j < singularvalcut[i].size(); j++) {
      Si += pow(abs(singularvalcut[i](j)), 2) * log(pow(abs(singularvalcut[i](j)), 2));

    }
    entropy(i) = (-1)*Si;
    Si=0;
  }

  return entropy;
  }

template <typename T>
double truncerr( Vector<T, Dynamic> sinvalj, Vector<T, Dynamic> sinvalcutj){
  int mj = sinvalcutj.size();
  double truncerr = 0;
  for (size_t i=mj; i < sinvalj.size(); i++) {
    truncerr += pow(abs(sinvalj(i)),2);
  }

  return truncerr;
}


template <typename T>
T coefofstate(std::vector<Matrix<T, Dynamic, Dynamic>> &matrices, VectorXi &state,int Nspin, int d){

  Matrix<T, Dynamic, Dynamic> A(matrices[0].rows(), matrices[0].cols());
  A = matrices[0](seq(state(0),last,d),all);
  Matrix<T, Dynamic, Dynamic> matsmall;

  for (size_t j = 1; j < (Nspin-1); j++) {
    matsmall = matrices[j](seq(state(j),last,d),all);
    Matrix<T, Dynamic, Dynamic> temp = A*matsmall;
    A.resize(temp.rows(), temp.cols());
    A = temp;

  }

  matsmall = matrices[Nspin-1](seq(state(Nspin-1),last,d),all);

  return (A*matsmall)(0,0);

}

template<typename T, int Nspin >
Vector<T, Dynamic> veckoeffrmps(std::vector<Matrix<T, Dynamic, Dynamic>> &matrices, int d){
  Vector<T, Dynamic> koef((int) pow(d, Nspin));

  for (size_t i = 0; i < pow(d, Nspin); i++) {
    auto bind = bitset<Nspin>(i);
    std::string str = bind.to_string();
    VectorXi state(Nspin);
    for (size_t j = 0; j < Nspin; j++) {
      state(j) = int(str[j])-48;
    }
    koef(i) = coefofstate<T>(matrices,state,Nspin, d);
  }

  return koef;
}

template <typename T>
std::vector<std::vector<Matrix<T, Dynamic, Dynamic>>>
MPSwB(Vector<T, Dynamic>  &stateref,int Nspin,
  int d, std::string odrezi, double treshold){
    auto Adecomp = MPSdecomp(stateref,Nspin,d, odrezi, treshold);

    for (size_t i = 1; i < Adecomp[0].size()-1; i++) {
      Vector<T, Dynamic> inv = 1.0/Adecomp[3][i-1].array();
      Matrix<T, Dynamic, Dynamic> diag = inv.asDiagonal();
      Adecomp[0][i](seq(0, last, 2), all ) = diag* Adecomp[0][i](seq(0, last, 2), all ).eval();
      Adecomp[0][i](seq(1, last, 2), all ) = diag* Adecomp[0][i](seq(1, last, 2), all ).eval();
    }

    Vector<T, Dynamic> inv = 1.0/Adecomp[3][Adecomp[0].size()-2].array();
    Adecomp[0][Adecomp[0].size()-1](0, all ) = inv.asDiagonal() *
     Adecomp[0][Adecomp[0].size()-1](0, all ).transpose().eval();
    Adecomp[0][Adecomp[0].size()-1](1, all ) = inv.asDiagonal() *
    Adecomp[0][Adecomp[0].size()-1](1, all ).transpose().eval();


    return {Adecomp[0],Adecomp[3]};
  }

template<typename T>
T coefofstateB(std::vector<Matrix<T, Dynamic, Dynamic>> &matrices,VectorXi &state,
  std::vector<Matrix<T, Dynamic, Dynamic>> & singularval, int Nspin,int d){
    Matrix<T, Dynamic, Dynamic> B(matrices[0].rows(), matrices[0].cols());
    B = matrices[0](seq(state(0),last,d),all);
    Matrix<T, Dynamic, Dynamic> matsmall;

    for (size_t j = 1; j < (Nspin-1); j++) {

      Matrix<T, Dynamic, Dynamic> diag = singularval[j-1].asDiagonal();
      matsmall = diag*matrices[j](seq(state(j),last,d),all).eval();
      B = B*matsmall;
      //B.resize(temp.rows(), temp.cols());
      //B = temp;

    }
    Matrix<T, Dynamic, Dynamic> diag = singularval[Nspin-2].asDiagonal();
    matsmall = diag*(matrices[Nspin-1](seq(state(Nspin-1),last,d), all ).transpose().eval());
    B = B*matsmall;
    return B(0,0);
  }

template<typename T, int Nspin>
Vector<T, Dynamic> veckoeffrmpsB(std::vector<Matrix<T, Dynamic, Dynamic>> &matrices,
  std::vector<Matrix<T, Dynamic, Dynamic>> & singularval, int d){
  Vector<T, Dynamic> koef((int) pow(d, Nspin));

  for (size_t i = 0; i < pow(d, Nspin); i++) {
    auto bind = bitset<Nspin>(i);
    std::string str = bind.to_string();
    VectorXi state(Nspin);
    for (size_t j = 0; j < Nspin; j++) {
      state(j) = int(str[j])-48;
    }
    koef(i) = coefofstateB<T>(matrices,state,singularval, Nspin, d);
  }

  return koef;
}

template<typename T>
Matrix<T, Dynamic, Dynamic> TijA(Matrix<T, Dynamic, Dynamic> &A1ij,
Matrix<T, Dynamic, Dynamic> &A2ij, int d){

  int rows = A1ij.rows()*A2ij.rows()/(pow(d,2));
  int cols = A1ij.cols()*A2ij.cols();
  Matrix<T, Dynamic, Dynamic> temp = MatrixXd::Zero(rows, cols);
  for (size_t j = 0; j < d; j++) {
    temp += kroneckerProduct(A1ij(seq(j, last, d), all).conjugate(), A2ij(seq(j, last, d), all));
  }

  return temp;
}

// THESE SCALAR PRODUCTS ARE MEMORY UNEFFICIENT!!

template<typename T>
double TscalprodA(std::vector<Matrix<T, Dynamic, Dynamic>> &A1ij,
std::vector<Matrix<T, Dynamic, Dynamic>> &A2ij, int d){

  auto prod = TijA<T>(A1ij[0], A2ij[0], 2);
  for (size_t i = 1; i < A1ij.size()-1; i++) {
    prod *= TijA<T>(A1ij[i], A2ij[i], 2);
  }
  //last
  prod *= TijA<T>(A1ij[A1ij.size()-1], A2ij[A1ij.size()-1], 2).transpose();

  return prod(0,0).real();
}

template<typename T>
Matrix<T, Dynamic, Dynamic> TijB(Matrix<T, Dynamic, Dynamic> &B1ij,
Matrix<T, Dynamic, Dynamic> &B2ij, Matrix<T, Dynamic, Dynamic> &lambda1,
Matrix<T, Dynamic, Dynamic> &lambda2 ,int d, std::string prva){

  std::vector<Matrix<T, Dynamic, Dynamic>> Tij(B1ij.size()-2);

  int rows = B1ij.rows()*B2ij.rows()/(pow(d,2));
  int cols = B1ij.cols()*B2ij.cols();
  Matrix<T, Dynamic, Dynamic> temp = MatrixXd::Zero(rows, cols);
  if (prva == "prva"){
    for (size_t j = 0; j < d; j++) {
      temp += kroneckerProduct(B1ij(seq(j, last, d), all).conjugate(),
      B2ij(seq(j, last, d), all));
    }
  }
  else if (prva == "zadnja"){
    for (size_t j = 0; j < d; j++) {
      temp += kroneckerProduct((lambda1.asDiagonal()*(B1ij(seq(j, last, d), all).transpose())).conjugate(),
      lambda2.asDiagonal()*(B2ij(seq(j, last, d), all).transpose()));
    }
  }

  else{
    for (size_t j = 0; j < d; j++) {
      temp += kroneckerProduct((lambda1.asDiagonal()*B1ij(seq(j, last, d), all)).conjugate(),
      lambda2.asDiagonal()*B2ij(seq(j, last, d), all));
    }
  }


  return temp;
}



template<typename T>
double TscalprodB(std::vector<Matrix<T, Dynamic, Dynamic>> &B1ij,
  std::vector<Matrix<T, Dynamic, Dynamic>> &B2ij,
  std::vector<Matrix<T, Dynamic, Dynamic>> &lambda1ij,
  std::vector<Matrix<T, Dynamic, Dynamic>> &lambda2ij,int d){

  auto prod = TijB<T>(B1ij[0], B2ij[0],lambda1ij[0], lambda2ij[0], 2, "prva");
  for (size_t i = 1; i < B1ij.size()-1; i++) {
    prod *= TijB<T>(B1ij[i], B2ij[i],lambda1ij[i-1], lambda2ij[i-1], 2, "splosno");
  }
  //last:
  int i = B1ij.size()-1;
  prod *= TijB<T>(B1ij[i], B2ij[i],lambda1ij[i-1], lambda2ij[i-1], 2, "zadnja").transpose();

  return prod(0,0).real();
}

template<typename T>
Matrix<T, Dynamic, Dynamic> Vj(Matrix<T, Dynamic, Dynamic> &B1ij,
Matrix<T, Dynamic, Dynamic> &B2ij, Matrix<T, Dynamic, Dynamic> &lambda1,
Matrix<T, Dynamic, Dynamic> &lambda2 ,int d, std::string prva){

  std::vector<Matrix<T, Dynamic, Dynamic>> Tij(B1ij.size()-2);

  int rows = B1ij.rows()*B2ij.rows()/(pow(d,2));
  int cols = B1ij.cols()*B2ij.cols();
  Matrix<T, Dynamic, Dynamic> temp = MatrixXd::Zero(rows, cols);
  if (prva == "prva"){
      temp = kroneckerProduct(B1ij(seq(0, last, d), all).conjugate(),
      B2ij(seq(0, last, d), all))- kroneckerProduct(B1ij(seq(1, last, d), all).conjugate(),
      B2ij(seq(1, last, d), all));

  }

  else if (prva == "zadnja"){
    temp = kroneckerProduct((lambda1.asDiagonal()*(B1ij(seq(0, last, d), all).transpose())).conjugate(),
    lambda2.asDiagonal()*(B2ij(seq(0, last, d), all).transpose())) -
    kroneckerProduct((lambda1.asDiagonal()*(B1ij(seq(1, last, d), all).transpose())).conjugate(),
    lambda2.asDiagonal()*(B2ij(seq(1, last, d), all).transpose()));

  }


  else{
      temp = kroneckerProduct((lambda1.asDiagonal()*B1ij(seq(0, last, d), all)).conjugate(),
      lambda2.asDiagonal()*B2ij(seq(0, last, d), all)) -
      kroneckerProduct((lambda1.asDiagonal()*B1ij(seq(1, last, d), all)).conjugate(),
      lambda2.asDiagonal()*B2ij(seq(1, last, d), all));
  }



  return temp;
}

template<typename T>
double sigmaksigmaj(std::vector<Matrix<T, Dynamic, Dynamic>> &B1ij,
  std::vector<Matrix<T, Dynamic, Dynamic>> &B2ij,
  std::vector<Matrix<T, Dynamic, Dynamic>> &lambda1ij,
  std::vector<Matrix<T, Dynamic, Dynamic>> &lambda2ij,int d, int j, int k){

  Matrix<T, Dynamic, Dynamic> res;
  if((k==0 && j==0) || (j !=0 && k !=0)){
    res = TijB<T>(B1ij[0], B2ij[0],lambda1ij[0], lambda2ij[0], 2, "prva");
  }
  else {
    res = Vj<T>(B1ij[0], B2ij[0],lambda1ij[0], lambda2ij[0], 2, "prva");
  }
  for (size_t i = 1; i < B1ij.size()-1; i++) {
    if((k==i && j==i) || (j !=i && k !=i)){
      res *= TijB<T>(B1ij[i], B2ij[i],lambda1ij[i-1], lambda2ij[i-1], 2, "splosno");
    }
    else{
      res *= Vj<T>(B1ij[i], B2ij[i],lambda1ij[i-1], lambda2ij[i-1], 2, "splosno");
    }
  }
  //last
  int i = B1ij.size()-1;
  if((k==i && j==i) || (j !=i && k !=i)){
    res *= TijB<T>(B1ij[i], B2ij[i],lambda1ij[i-1], lambda2ij[i-1], 2, "zadnja");
  }
  else{
    res *= Vj<T>(B1ij[i], B2ij[i],lambda1ij[i-1], lambda2ij[i-1], 2, "zadnja");
  }
  return res(0,0).real();
}

/////////////////////////
//helper function:
//non-compact bipartition
/*
Vector2i ABABshuffleind(std::string & bit, int N){
  std::string levo = "";
  std::string desno = "";
  for (int i = 0; i <N ; i++) {
    if(i % 2 ==0){
      levo += bit[i];
    }
    else{
      desno += bit[i];
    }
  }
  Vector2i res;

  res(1) = stoi(levo,nullptr, 2);
  res(0) = stoi(desno,nullptr, 2);

  return res;
}

Vector2i AABBshuffleind(std::string & bit, int N){
  std::string levo = "";
  std::string desno = "";
  int left = 1;
  int stev = 0;
  for (int i = 0; i <N ; i++) {
    if(stev==2 && left==1){
      left=0;
      stev =0;
    }
    if(stev==2 && left==0){
      left=1;
      stev = 0;
    }
    stev +=1;

    if(left==1){
      levo += bit[i];

    }
    else{
      desno += bit[i];
    }
  }

  Vector2i res;

  res(1) = stoi(levo,nullptr, 2);
  res(0) = stoi(desno,nullptr, 2);

  return res;
}

template<typename T, int Nspin>
Matrix<T, Dynamic, Dynamic> shuffleAABBpart(Vector<T, Dynamic> &state){
  int Nhalf = pow(2, Nspin/2);
  Vector2i vecind;
  Matrix<T,Dynamic, Dynamic> res(Nhalf, Nhalf);
  if (Nspin==6){
    res.resize((int)pow(2, 2),(int)pow(2, 4));
  }
  else if (Nspin==10){
    res.resize((int)pow(2, 4),(int)pow(2, 6));
  }
  else if (Nspin==14){
    res.resize((int)pow(2, 6),(int)pow(2, 8));
  }

  for (unsigned int i = 0; i < state.size(); i++) {
    auto bind = bitset<Nspin>(i);
    std::string bitstr = bind.to_string();
    vecind = AABBshuffleind(bitstr, Nspin);
    res(vecind(0), vecind(1)) = state(i);

  }

  return res;
}
template<typename T, int Nspin>
Matrix<T, Dynamic, Dynamic> shuffleABABpart(Vector<T, Dynamic> &state){
  constexpr int Nhalf = pow(2, Nspin/2);
  Matrix<T,Nhalf, Nhalf> res;

  for (unsigned int i = 0; i < state.size(); i++) {
    auto bind = bitset<Nspin>(i);
    std::string bitstr = bind.to_string();
    Vector2i vecind= ABABshuffleind(bitstr, Nspin);
    res(vecind(0), vecind(1)) = state(i);

  }
  return res;
}

template<typename T, int Nspin>
double ABABentropysim(Vector<T, Dynamic> &state, double treshold){
  //Matrix<T, Dynamic, Dynamic> mat = shuffleABABpart<T, Nspin>(state);
  Matrix<T, Dynamic, Dynamic> mat = shuffleAABBpart<T, Nspin>(state);
  BDCSVD<Matrix<T, Dynamic, Dynamic>> bdcsvd( mat,  ComputeThinU | ComputeThinV);
  bdcsvd.setThreshold(treshold);
  int mi = bdcsvd.rank();
  double S=0;
  for (size_t i = 0; i < mi; i++) {
    if (not(isnan(bdcsvd.singularValues()(i)))){
      S += pow(abs(bdcsvd.singularValues()(i)), 2) * log(pow(abs(bdcsvd.singularValues()(i)), 2));
    }
  }
  return -1.0*S;
}
*/










#endif
