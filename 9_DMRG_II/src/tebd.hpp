#ifndef TEBD_HPP
#define TEBD_HPP
//#include "../8_DMRG_I/dmrg.hpp"
#include "../8_DMRG_I/mps_algorithm.hpp"


template<typename T>
Matrix<T, Dynamic, Dynamic>  U(T z, double a){
  Matrix<T, Dynamic,Dynamic> Umat = Matrix<T, Dynamic, Dynamic>::Zero(4,4);
  Umat(0,0) =Umat(3,3) =  exp(z*a);
  Umat(1,1) = Umat(2,2) = exp(-z*a)*cosh(2.0*z*a);
  Umat(1,2) = Umat(2,1) = exp(-z*a)*sinh(2.0*z*a);

  return Umat;
}

template<typename T>
void localTEBDstep(std::vector<Matrix<T, Dynamic, Dynamic>> & Bij,
std::vector<Matrix<T, Dynamic, Dynamic>> &lambde, Matrix<T, Dynamic, Dynamic> &U ,int j, double treshold,
double &localerror) {

  Matrix<T, Dynamic, Dynamic> B00 = U(0,0)* Bij[j](seq(0, last, 2), all).eval()*lambde[j].asDiagonal()* Bij[j+1](seq(0, last, 2), all).eval();
  Matrix<T, Dynamic, Dynamic> B01 = U(1,1)* Bij[j](seq(0, last, 2), all).eval()*lambde[j].asDiagonal()* Bij[j+1](seq(1, last, 2), all).eval() +
  U(1,2)* Bij[j](seq(1, last, 2), all).eval()*lambde[j].asDiagonal()* Bij[j+1](seq(0, last, 2), all).eval();
  Matrix<T, Dynamic, Dynamic> B10  = U(2,2)* Bij[j](seq(1, last, 2), all).eval()*lambde[j].asDiagonal()* Bij[j+1](seq(0, last, 2), all).eval() +
  U(2,1)* Bij[j](seq(0, last, 2), all).eval()*lambde[j].asDiagonal()* Bij[j+1](seq(1, last, 2), all).eval();
  Matrix<T, Dynamic, Dynamic> B11 = U(3,3)* Bij[j](seq(1, last, 2), all).eval()*lambde[j].asDiagonal()* Bij[j+1](seq(1, last, 2), all).eval();



  Matrix<T, Dynamic, Dynamic> Q;
  Matrix<T, Dynamic, Dynamic> Q0;
  Matrix<T, Dynamic, Dynamic> Q1;

  if (j==0){
    Q = Matrix<T, Dynamic, Dynamic>::Zero(2,(lambde[j+1].size()*2));
    Q0 = Matrix<T, Dynamic, Dynamic>::Zero(2,(lambde[j+1].size()));
    Q1 = Matrix<T, Dynamic, Dynamic>::Zero(2,(lambde[j+1].size()));


    Q0(seq(0, last, 2), all) = B00*lambde[j+1].asDiagonal();
    Q0(seq(1, last, 2), all) = B10*lambde[j+1].asDiagonal();
    Q1(seq(0, last, 2), all) = B01*lambde[j+1].asDiagonal();
    Q1(seq(1, last, 2), all) = B11*lambde[j+1].asDiagonal();
  }

  else if (j==(Bij.size()-2)){


    Q = Matrix<T, Dynamic, Dynamic>::Zero((lambde[j-1].size()*2),2);
    Q0 = Matrix<T, Dynamic, Dynamic>::Zero((lambde[j-1].size())*2,1);
    Q1 = Matrix<T, Dynamic, Dynamic>::Zero((lambde[j-1].size())*2,1);


    Q0(seq(0, last, 2), all) = lambde[j-1].asDiagonal()*B00;
    Q0(seq(1, last, 2), all) = lambde[j-1].asDiagonal()*B10;
    Q1(seq(0, last, 2), all) = lambde[j-1].asDiagonal()*B01;
    Q1(seq(1, last, 2), all) = lambde[j-1].asDiagonal()*B11;

  }
  else{

    Q = Matrix<T, Dynamic, Dynamic>::Zero((lambde[j-1].size()*2),(lambde[j+1].size()*2));
    Q0 = Matrix<T, Dynamic, Dynamic>::Zero((lambde[j-1].size()*2),(lambde[j+1].size()));
    Q1 = Matrix<T, Dynamic, Dynamic>::Zero((lambde[j-1].size()*2),(lambde[j+1].size()));

    Q0(seq(0, last, 2), all) = (lambde[j-1].asDiagonal()*B00)*lambde[j+1].asDiagonal();
    Q0(seq(1, last, 2), all) = (lambde[j-1].asDiagonal()*B10)*lambde[j+1].asDiagonal();
    Q1(seq(0, last, 2), all) = (lambde[j-1].asDiagonal()*B01)*lambde[j+1].asDiagonal();
    Q1(seq(1, last, 2), all) = (lambde[j-1].asDiagonal()*B11)*lambde[j+1].asDiagonal();
  }



  Q(all, seq(0, last,2)) = Q0;
  Q(all, seq(1, last,2)) = Q1;




  BDCSVD<Matrix<T, Dynamic, Dynamic>> bdcsvd1( Q,  ComputeThinU | ComputeThinV);


  int M = lambde[j].size();
  double maxsing = bdcsvd1.singularValues().maxCoeff();
  bdcsvd1.setThreshold(treshold/abs(maxsing));
  int mi = bdcsvd1.singularValues().size();

  int mi2 = bdcsvd1.rank();
  if (mi2 >= M){
    mi = M;
  }

  else{
    std::cout << "manjsi, tezave!" << '\n';
    mi = mi2;
  }


  if (bdcsvd1.singularValues()(0) < pow(10, -8)){
    //std::cout <<"majhne lastne vrednosti"<< bdcsvd1.singularValues() << '\n';
    //exit(1);
  }

  for (size_t i =(mi+1); i < (bdcsvd1.singularValues().size()); i++) {
    localerror += pow(abs(bdcsvd1.singularValues()(i)),2);
  }


  Vector<T, Dynamic> invsingjminus1;
  Vector<T, Dynamic> invsingjplus1;
  Matrix<T, Dynamic, Dynamic> Vdagger = bdcsvd1.matrixV().adjoint();
  Matrix<T, Dynamic, Dynamic> Bjnew = Matrix<T, Dynamic, Dynamic>::Zero(bdcsvd1.matrixU().rows()
  ,bdcsvd1.matrixU()(all, seq(0, mi-1)).cols());
  Matrix<T, Dynamic, Dynamic> Bjplus1new = Matrix<T, Dynamic, Dynamic>::Zero(
    (mi)*2, Vdagger.cols()/2);

  if (j==0){
    invsingjplus1  = 1.0/(lambde[j+1].array());
    Matrix<T, Dynamic, Dynamic> invplsdiag = invsingjplus1.asDiagonal();
    Bjnew= bdcsvd1.matrixU()(all, seq(0, mi-1));
    Bjplus1new(seq(0, last,2),all)= Vdagger(seq(0,mi-1), seq(0, last,2)).eval()
    * invplsdiag;
    Bjplus1new(seq(1, last,2), all )= Vdagger(seq(0,mi-1), seq(1, last,2)).eval()
    * invplsdiag;
  }


  else if(j==(Bij.size()-2)){

    invsingjminus1 = 1.0/(lambde[j-1].array());
    Matrix<T, Dynamic, Dynamic> invminusdiag = invsingjminus1.asDiagonal();
    Bjnew(seq(0,last, 2), all)= invminusdiag*(bdcsvd1.matrixU()(seq(0,last, 2), seq(0,mi-1))).eval();
    Bjnew(seq(1,last, 2), all)= invminusdiag*(bdcsvd1.matrixU()(seq(1,last, 2), seq(0,mi-1))).eval();
    Bjplus1new = Matrix<T, Dynamic, Dynamic>::Zero(
      (mi), Vdagger.cols());
    Bjplus1new(seq(0,last, 2), all)= Vdagger(all, 0).transpose();
    Bjplus1new(seq(1,last, 2), all)= Vdagger(all, 1).transpose();
  }
  else {
    invsingjminus1  = 1.0/(lambde[j-1].array());
    invsingjplus1  = 1.0/(lambde[j+1].array());
    Matrix<T, Dynamic, Dynamic> invminusdiag = invsingjminus1.asDiagonal();
    Matrix<T, Dynamic, Dynamic> invplsdiag = invsingjplus1.asDiagonal();

    Bjnew(seq(0,last, 2),all)= invminusdiag*((bdcsvd1.matrixU()(seq(0,last, 2), seq(0,mi-1)))).eval();
    Bjnew(seq(1,last, 2),all)= invminusdiag*((bdcsvd1.matrixU()(seq(1,last, 2), seq(0,mi-1)))).eval();


    Bjplus1new(seq(0,last, 2),all)= Vdagger(seq(0,mi-1),seq(0,last, 2)).eval()*invplsdiag;
    Bjplus1new(seq(1,last, 2),all)= Vdagger(seq(0,mi-1),seq(1,last, 2)).eval()* invplsdiag;
  }

  VectorXd singvalj = bdcsvd1.singularValues()(seq(0,mi-1));

  Bij[j] = Bjnew;
  Bij[j+1] = Bjplus1new;
  lambde[j] = singvalj;

}

template<typename T>
void U_even(std::vector<Matrix<T, Dynamic, Dynamic>> & Bij,
std::vector<Matrix<T, Dynamic, Dynamic>> &lambde, Matrix<T, Dynamic, Dynamic> &U
 , double treshold, double &localerror, int Nspin){
  for (size_t i = 0; i < Nspin/2; i++) {
    int i_ev = 2*i;
    localTEBDstep<T>(Bij, lambde, U, i_ev, treshold, localerror);
  }
}

template<typename T>
void U_odd(std::vector<Matrix<T, Dynamic, Dynamic>> & Bij,
std::vector<Matrix<T, Dynamic, Dynamic>> &lambde, Matrix<T, Dynamic, Dynamic> &U
 , double treshold, double &localerror, int Nspin){
  for (size_t i = 0; i < Nspin/2-1; i++) {
    int i_odd = 2*i+1;
    localTEBDstep<T>(Bij, lambde, U, i_odd, treshold, localerror);

  }
}


template<typename T>
void S2(std::vector<Matrix<T, Dynamic, Dynamic>> & Bij,
std::vector<Matrix<T, Dynamic, Dynamic>> &lambde
 , double treshold, double &localerror, complex<double> z , int Nspin){

  auto U1 = U<complex<double>>(z, 0.5);
  auto U2 = U<complex<double>>(z, 1.0);

  U_even<T>(Bij, lambde, U1, treshold, localerror, Nspin);
  U_odd<T>(Bij, lambde, U2, treshold, localerror, Nspin);
  U_even<T>(Bij, lambde, U1, treshold, localerror, Nspin);

}

template<typename T>
void S4(std::vector<Matrix<T, Dynamic, Dynamic>> & Bij,
std::vector<Matrix<T, Dynamic, Dynamic>> &lambde
 , double treshold, double &localerror, complex<double> z, int Nspin){

  double k1 = 1.0/(2-pow(2.0,(1.0/3.0)));
  double k0 = - pow(2.0,(1.0/3.0))*k1;

  S2<T>(Bij,lambde,treshold,localerror, k1*z, Nspin);
  S2<T>(Bij,lambde,treshold,localerror, k0*z, Nspin);
  S2<T>(Bij,lambde,treshold,localerror, k1*z, Nspin);

 }

template<typename T>
void suz_trotter(std::vector<Matrix<T, Dynamic, Dynamic>> & Bij,
std::vector<Matrix<T, Dynamic, Dynamic>> &lambde
 , double treshold, double &localerror, complex<double> z, int Nspin, int koraki,
  std::string metoda) {

  complex<double> curz = 0.0;

    for (size_t i = 0; i < koraki; i++) {
      curz +=abs(z.real());
      if (metoda =="S2"){
        S2<T>(Bij,lambde,treshold,localerror, z, Nspin);
        //std::cout << lambde[10] << '\n';
      }
      else if (metoda == "S4"){
        S4<T>(Bij,lambde,treshold,localerror, z, Nspin);
      }
      else {
        std::cout << "Izberi pravo metodo" << '\n';
        exit(1);
      }
    }
}


template<typename T>
std::vector<std::vector<Matrix<T, Dynamic, Dynamic>>>
 MPSrandomstate(int Nspin, int cutoffsize, std::string normalize){

  if (int(log2(cutoffsize)) >= Nspin/2){
    cutoffsize = pow(2, Nspin/2);
  }

  VectorXi dimrows(Nspin);
  VectorXi dimcols(Nspin);
  for (size_t i = 0; i < int(log2(cutoffsize)); i++) {
    dimcols(i) = pow(2, i+1);
  }
  int Nsame = Nspin - 2*int(log2(cutoffsize))-1;


  for (size_t j = int(log2(cutoffsize)); j < (int(log2(cutoffsize))+Nsame); j++) {
    dimcols(j) = cutoffsize;
  }
  dimcols(seq((int(log2(cutoffsize))+Nsame), Nspin-2)) = dimcols(seq(0,int(log2(cutoffsize))-1)).reverse();
  dimcols(Nspin-1) = 2;
  dimrows(0) = 1;
  for (size_t j = 1; j < Nspin; j++){
    dimrows(j) = dimcols(j-1);
  }

  dimrows(seq(0, last-1)) *=2;


  std::vector<Matrix<T, Dynamic, Dynamic>> Aij(Nspin);
  std::vector<Matrix<T, Dynamic, Dynamic>> lambdeId(Nspin-1);
  for (size_t j = 0; j < Nspin; j++) {
    Aij[j] = Matrix<T, Dynamic, Dynamic>::Random(dimrows(j), dimcols(j));

    if (normalize=="True"){
      Aij[j].normalize();
    }

    if (j< Nspin-1){
      lambdeId[j] = MatrixXd::Constant( dimcols(j),1, 1.0);
    }

  }

  return {Aij, lambdeId};

}


template<typename T>
T collapse_state(std::vector<Matrix<T, Dynamic, Dynamic>> & Bij,
std::vector<Matrix<T, Dynamic, Dynamic>> &lambde, VectorXi &randstate){
  MatrixXcd mat = Bij[0](seq(randstate(0),last,2), all);
  for (size_t i = 1; i < randstate.size()-1; i++) {
    mat = mat*lambde[i-1].asDiagonal()*Bij[i](seq(randstate(i),last,2), all).eval();
  }
  int i = randstate.size()-1;
  mat = mat*lambde[i-1].asDiagonal()*((Bij[i](seq(randstate(i),last,2), all).eval()).transpose());

  return mat(0,0);

}

template<typename T>
void sigma_zj_MPA
(std::vector<Matrix<T, Dynamic, Dynamic>> & Bij, int j){
  Bij[j](seq(1,last,2), all) *= (-1.0);
}

template<typename T>
std::vector<std::vector<Matrix<T, Dynamic, Dynamic>>>
MPSdomainstate(int Nspin, int cutoffsize, std::string normalize, VectorXi domainconf, double ksi){

  if (Nspin % 2 == 1){
    std::cout << "Vstavi sodo število spinov" << '\n';
    exit(1);
  }
  if (domainconf.size() > Nspin/2 || domainconf.sum() > Nspin || Nspin % domainconf.size() != 0 ){
    std::cout << "Konfiguracija domen ni prava" << '\n';
    exit(1);
  }

  if (int(log2(cutoffsize)) >= Nspin/2){
    cutoffsize = pow(2, Nspin/2);
  }

  VectorXi dimrows(Nspin);
  VectorXi dimcols(Nspin);
  for (size_t i = 0; i < int(log2(cutoffsize)); i++) {
    dimcols(i) = pow(2, i+1);
  }
  int Nsame = Nspin - 2*int(log2(cutoffsize))-1;


  for (size_t j = int(log2(cutoffsize)); j < (int(log2(cutoffsize))+Nsame); j++) {
    dimcols(j) = cutoffsize;
  }
  dimcols(seq((int(log2(cutoffsize))+Nsame), Nspin-2)) = dimcols(seq(0,int(log2(cutoffsize))-1)).reverse();
  dimcols(Nspin-1) = 2;
  dimrows(0) = 1;
  for (size_t j = 1; j < Nspin; j++){
    dimrows(j) = dimcols(j-1);
  }

  dimrows(seq(0, last-1)) *=2;


  std::vector<Matrix<T, Dynamic, Dynamic>> Aij(Nspin);
  std::vector<Matrix<T, Dynamic, Dynamic>> lambdeId(Nspin-1);

  int spindir = 1;
  int currind = 0;
  for (size_t k = 0; k < domainconf.size(); k++) {
    for (size_t i = 0; i < domainconf(k); i++) {
      Array<T, Dynamic, Dynamic> temp = Array<T, Dynamic, Dynamic>::Random(dimrows(currind), dimcols(currind));
      temp(seq(spindir,last, 2), all) *= ksi;
      Aij[currind] = temp.matrix();

      if (normalize=="True"){
        Aij[currind].normalize();
      }

      if (currind < Nspin-1){
        lambdeId[currind] = MatrixXd::Constant( dimcols(currind),1, 1.0);
        if (currind >=1) {
          MatrixXcd id = lambdeId[currind-1].asDiagonal();

        }
      }

      currind+=1;
    }
    if (spindir == 1){
      spindir = 0;
    }
    else {
      spindir = 1;
    }

  }



  return {Aij, lambdeId};
}






























#endif
