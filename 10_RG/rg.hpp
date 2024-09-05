#ifndef RG_HPP
#define RG_HPP
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/CXX11/Tensor>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <cmath>


using namespace std;
using namespace Eigen;

namespace squarelatt
{
  int qDitToDecimal(int q, const std::vector<int>& digits) {
      int decimalValue = 0;
      int power = digits.size() - 1;

      for (int digit : digits) {
          decimalValue += digit * pow(q, power);
          power--;
      }

      return decimalValue;
  }

  std::vector<int> twodecToBaseQ(int decimalValue, int q) {
      std::vector<int> digits;
      digits.push_back(decimalValue/q);
      digits.push_back(decimalValue%q);

      return digits;
  }


  double deltaf(int s1, int s2){
    if(s1==s2){
      return 1.0;
    }
    else {
      return 0.0;
    }
  }

  Tensor<double, 4> A0tensor(int q, double beta, double J){
    Tensor<double, 4> A(q,q,q,q);
    for (int l = 0; l < q; l++) {
      for (int r = 0; r < q; r++) {
        for (int t = 0; t < q; t++) {
          for (int d = 0; d < q; d++) {
            A(l,r,t,d) = exp((beta*J/2.0)*
            (deltaf(t+1, r+1) + deltaf(l+1, t+1) + deltaf(d+1, l+1) + deltaf(d+1, r+1)));

          }
        }
      }
    }

    return (A);

  }
  /*
  vector<vector<Tensor<double, 4>>> make_lattice(int N, Tensor<double, 4> ten){
    int s = int(pow(2, N));
    vector<vector<Tensor<double, 4>>> mat(s, vector<Tensor<double, 4>>(s));
    for (int i = 0; i < s; i++) {
      for (int j = 0; j < s; j++) {
        mat[i][j] = ten;
      }
    }
    return mat;

  }
  */

  Tensor<double, 3> Convmat(MatrixXd &mat, std::string F){
      int qnew = 0;
      int qold = 0;
      Tensor<double, 3> T;
      if (F=="F12"){
        qnew = mat.cols();
        qold = int(sqrt(mat.rows()));
        T.resize(qold, qold, qnew);
        for (int i = 0; i < qnew; i++) {
          for (int k = 0; k < int(pow(qold, 2)); k++) {
            std::vector<int> vec = twodecToBaseQ(k, qold);
            int prvi = vec[0];
            int drugi = vec[1];
            T(prvi, drugi, i) = mat(k, i);
          }
        }
      }

      else if (F=="F34"){
        qnew = mat.rows();
        qold = int(sqrt(mat.cols()));
        T.resize(qnew, qold, qold);

        for (int i = 0; i < qnew; i++) {
          for (int k = 0; k < int(pow(qold, 2)); k++) {
            std::vector<int> vec = twodecToBaseQ(k, qold);
            int prvi = vec[0];
            int drugi = vec[1];
            T(i,prvi, drugi) = mat(i, k);
          }
        }

      }


      else {
        std::cout << "vpisi pravi parameter" << '\n';
        exit(1);
      }

      return T;

  }



  MatrixXd AtoFF( Tensor<double, 4> &ten, std::string type){
    int dim = ten.dimensions()[0];
    MatrixXd mat(int(pow(dim,2)), int(pow(dim, 2)));
    int row=0;
    int col = 0;
    if (type=="sz"){
      for (int l = 0; l < dim; l++) {
        for (int d = 0; d < dim; d++) {
          for (int r = 0; r < dim; r++) {
            for (int t = 0; t < dim; t++) {
              mat(qDitToDecimal(dim, {l,t}), qDitToDecimal(dim,{d,r})) = ten(l,r,t,d);
            }
          }
        }
      }
    }
    else if(type=="sv"){
      for (int l = 0; l < dim; l++) {
        for (int t = 0; t < dim; t++) {
          for (int r = 0; r < dim; r++) {
            for (int d = 0; d < dim; d++) {
              mat(qDitToDecimal(dim, {t,r}), qDitToDecimal(dim,{l,d})) = ten(l,r,t,d);
            }
          }
        }
      }
    }
    else {
      std::cout << "vpisi pravi razcep!" << '\n';
      exit(1);
    }
    return mat;

  }

  std::vector<MatrixXd> computeFF(MatrixXd &FF, int maxdim){


    BDCSVD<MatrixXd> bdcsvd( FF,  ComputeThinU | ComputeThinV);
    int info = bdcsvd.info();
    /*
    // Check the status
    if (info == Eigen::Success) {
      // The SVD method terminated successfully.
      std::cout << "sucess" << '\n';
    } else if (info == Eigen::NumericalIssue) {
      std::cout << "numerical issue" << '\n';
      // The SVD method failed due to a numerical issue.
    } else if (info == Eigen::InvalidInput) {
      std::cout << "invalid input" << '\n';
    } else if (info == Eigen::NoConvergence) {
      std::cout << "did not converge" << '\n';
    }
    */

    double maxsing = bdcsvd.singularValues().maxCoeff();

    //bdcsvd.setThreshold(pow(10, -10)/abs(maxsing));
    //int mi = bdcsvd.rank();
    int mi = bdcsvd.singularValues().size();
    if (mi > maxdim){
      mi = maxdim;
    }

    VectorXd sqrsingvec = bdcsvd.singularValues().cwiseSqrt();
    sqrsingvec = sqrsingvec(seq(0, mi-1)).eval();
    MatrixXd F34 = bdcsvd.matrixV().transpose()(seq(0, mi-1), all).eval();
    MatrixXd F12 = bdcsvd.matrixU()(all,seq(0, mi-1)).eval();
    F34 = (sqrsingvec.asDiagonal() * F34).eval();

    F12 = (F12*sqrsingvec.asDiagonal()).eval();

    std::vector<MatrixXd> stor = {F12, F34};

    return stor;
  }

  double Anewcomp(int lnew,int rnew,int tnew, int dnew, MatrixXd &F3, MatrixXd &F4,
  MatrixXd &F2, MatrixXd &F1){


    int size = int(sqrt(F3.cols()));

    double sum = 0;
    for (int x = 0; x < size; x++) {
      for (int y = 0; y < size; y++) {
        for (int z = 0; z < size; z++) {
          for (int w = 0; w < size; w++) {
            sum += F3(lnew, qDitToDecimal(size, {z,w})) * F4(tnew,qDitToDecimal(size, {w,y}))*
            F1(qDitToDecimal(size, {x,y}),rnew)*F2(qDitToDecimal(size, {z,x}),dnew);

          }
        }
      }
    }
  return sum;
  }

  Tensor<double, 4> Anew(MatrixXd &F3, MatrixXd &F4,MatrixXd &F2, MatrixXd &F1){
    int dim = F1.cols();
    int oldim = F1.rows();



    Tensor<double, 3> TF1 = Convmat(F1, "F12");
    Tensor<double, 3> TF2 = Convmat(F2, "F12");
    Tensor<double, 3> TF3 = Convmat(F3, "F34");
    Tensor<double, 3> TF4 = Convmat(F4, "F34");


    //TENSOR MAP DOES NOT GIVE RIGHT RESULT, I THINK DUE TO COL-MAJOR ORDER
    /*
    Eigen::array<int, 3> f12info{{int(sqrt(F1.rows())),int(sqrt(F1.rows())),int(F1.cols())}};
    Eigen::array<int, 3> f34info{{int(F3.rows()),int(sqrt(F3.cols())),int(sqrt(F3.cols()))}};

    Tensor<double, 3> TF1= TensorMap<Eigen::Tensor<double, 3>>(F1.data(), f12info);
    Tensor<double, 3> TF2= TensorMap<Eigen::Tensor<double, 3>>(F2.data(), f12info);
    Tensor<double, 3> TF3= TensorMap<Eigen::Tensor<double, 3>>(F3.data(), f34info);
    Tensor<double, 3> TF4= TensorMap<Eigen::Tensor<double, 3>>(F4.data(), f34info);
    */

    Eigen::array<IndexPair<int>, 1> fst = { IndexPair<int>(1, 2) };
    Tensor<double, 4> Temp = TF1.contract(TF4, fst);
    Eigen::array<IndexPair<int>, 1> scnd = { IndexPair<int>(2, 3) };
    Tensor<double, 5> Temp2 = TF3.contract(Temp, scnd);
    Eigen::array<IndexPair<int>, 2> last = { IndexPair<int>(1, 0),IndexPair<int>(2, 1) };
    Temp = Temp2.contract(TF2, last).eval();

    auto arr = Temp.dimensions();
    std::cout << "[" <<arr[0]<< ","<<arr[1]<<","<<arr[2]<<","<<arr[3]<<"]" <<'\n';

    /*
    Tensor<double, 4> An(dim,dim,dim,dim);
    for (int l = 0; l < dim; l++) {
      for (int r = 0; r < dim; r++) {
        for (int t = 0; t < dim; t++) {
          for (int d = 0; d < dim; d++) {
            An(l,r,t,d) = Anewcomp(l,r,t,d, F3,F4,F2,F1);
          }
        }
      }
    }
    auto arr = An.dimensions();
    std::cout << "[" <<arr[0]<< ","<<arr[1]<<","<<arr[2]<<","<<arr[3]<<"]" <<'\n';
    */
    return Temp;
  }

  double Ancontr(Tensor<double,4> & finaltens){
    int dim = finaltens.dimensions()[0];
    double sum = 0;
    for (int lr = 0; lr < dim; lr++) {
      for (int td = 0; td < dim; td++) {
        sum += finaltens(lr,lr,td,td);
      }
    }

    return sum;
  }


  Tensor<double, 4> RGstep(Tensor<double, 4> &Aj, int maxdim){
    MatrixXd matsz = AtoFF(Aj, "sz");
    MatrixXd matsv = AtoFF(Aj, "sv");



    std::vector<MatrixXd> vecsz = computeFF(matsz, maxdim);
    std::vector<MatrixXd> vecsv = computeFF(matsv, maxdim);


    return Anew(vecsz[1], vecsv[1], vecsv[0], vecsz[0]);

  }

  double calc_norm_factor(VectorXd &norms){
    double sum = 0;
    int N = norms.size();
    for (int i = 0; i < N; i++) {
      sum += pow(2, N-i) * log(norms(i));
    }
    return sum;
  }

  Vector2d Partitsum(int numq, double beta, double J, int steps, int maxdim){
    Tensor<double, 4> Aj = A0tensor(numq,beta,J);
    VectorXd norms(steps);
    int ind = 0;
    while ((steps) > pow(10, -4)) {
      Eigen::Map<Eigen::VectorXd> vec(Aj.data(), Aj.size());
      double norm = vec.norm();
      norms(ind) = norm;
      Aj = Aj/norm;
      Aj = RGstep(Aj, maxdim).eval();


      steps -= 1;
      ind +=1;
    }

    double res = Ancontr(Aj);
    double fact = calc_norm_factor(norms);
    Vector2d resall;
    resall << res, fact;
    return resall;
  }

  VectorXd lnZ_sweep_beta(int numq,double J, int steps, int maxdim){
      VectorXd betasm = VectorXd::LinSpaced(39, 0.05, 1.0);
      VectorXd betalg = VectorXd::LinSpaced(120, 1.025, 4.0);
      VectorXd beta(betasm.size()+ betalg.size());
      beta <<betasm , betalg;
      //std::cout << beta << '\n';
      VectorXd lnz(beta.size());
      for (size_t i = 0; i < beta.size(); i++) {
        //std::cout << (1.0*i/(beta.size()-1))*100 << '\n';
        auto resall = Partitsum(numq, beta(i), J, steps, maxdim);
        lnz(i) = log(resall(0))-resall(1);
      }
      return lnz;
    }


  void lnZ_sweep_beta_tofile(int numq,double J, int steps, int maxdim, std::string imedat){
      VectorXd betasm = VectorXd::LinSpaced(39, 0.05, 1.0);
      VectorXd betalg = VectorXd::LinSpaced(120, 1.025, 4.0);
      VectorXd beta(betasm.size()+ betalg.size());
      beta <<betasm , betalg;
      VectorXd lnz(beta.size());
      for (size_t i = 0; i < beta.size(); i++) {
        std::cout << (1.0*i/(beta.size()-1))*100 << '\n';
        auto resall = Partitsum(numq, beta(i), J, steps, maxdim);
        lnz(i) = log(resall(0))-resall(1);
      }

      std::ofstream outputFile(imedat, std::ios::out);
      outputFile << "beta, lnZ, db=" << beta(1)-beta(0) << "\n";
      for (size_t i = 0; i < beta.size(); i++) {
        outputFile << beta(i) << " " << lnz(i) << "\n";
      }

      outputFile.close();


    }
}

namespace honeycomb
{
  int qDitToDecimal(int q, const std::vector<int>& digits) {
      int decimalValue = 0;
      int power = digits.size() - 1;

      for (int digit : digits) {
          decimalValue += digit * pow(q, power);
          power--;
      }

      return decimalValue;
  }

  std::vector<int> twodecToBaseQ(int decimalValue, int q) {
      std::vector<int> digits;
      digits.push_back(decimalValue/q);
      digits.push_back(decimalValue%q);

      return digits;
  }


  double deltaf(int s1, int s2){
    if(s1==s2){
      return 1.0;
    }
    else {
      return 0.0;
    }
  }


  //INITIAL TENSORS T0tensor gives wrong result, must be other form
  Tensor<double, 3> T0tensor(int q, double beta, double J){
    Tensor<double, 3> T(q,q,q);
    for (int i = 0; i < q; i++) {
      for (int j = 0; j < q; j++) {
        for (int k = 0; k < q; k++) {
          T(i,j,k) = exp((beta*J/2.0)*
          (deltaf(i, j) + deltaf(j, k) + deltaf(k,i)));

        }
      }
    }
    return T;
  }


  std::vector<Tensor<double, 3>> T0ab(int q, double beta, double J){
    MatrixXd mat(q,q);
    for (size_t i = 0; i < q; i++) {
      for (size_t j = 0; j < q; j++) {
        mat(i,j) = exp((beta*J)*(deltaf(i,j)));
      }
    }

    BDCSVD<MatrixXd> bdcsvd( mat,  ComputeThinU | ComputeThinV);
    int mi = bdcsvd.singularValues().size();
    if (mi != q) {
      std::cout << "ERROR" << '\n';
      exit(1);
    }

    VectorXd sqrsingvec = bdcsvd.singularValues().cwiseSqrt();
    MatrixXd Qb = ( bdcsvd.matrixV()*sqrsingvec.asDiagonal()).eval();

    MatrixXd Qa = (bdcsvd.matrixU()*sqrsingvec.asDiagonal()).eval();

    Tensor<double, 3> Ta(mi,mi,mi);
    Tensor<double, 3> Tb(mi,mi,mi);
    for (int i = 0; i < mi; i++) {
      for (int j = 0; j < mi; j++) {
        for (int k = 0; k < mi; k++) {
          double tempsuma = 0;
          double tempsumb = 0;
          for (int l = 0; l < q; l++) {
            tempsuma += Qa(l,i)*Qa(l,j)*Qa(l,k);
            tempsumb += Qb(l,i)*Qb(l,j)*Qb(l,k);
          }
          Ta(i,j,k)=tempsuma;
          Tb(i,j,k)=tempsumb;
        }
      }
    }
    std::vector<Tensor<double,3>> vec = {Ta,Tb};
    return vec;
  }


  Tensor<double, 4> Aijkl(Tensor<double, 3> &Ta, Tensor<double,3> &Tb){
    Eigen::array<IndexPair<int>, 1> fst = { IndexPair<int>(2, 2) };
    Tensor<double, 4> A = Ta.contract(Tb, fst);

    return A;
  }

  MatrixXd AtoM(Tensor<double, 4> & Aijkl){
    int dim = Aijkl.dimensions()[0];

    MatrixXd mat(int(pow(dim,2)), int(pow(dim, 2)));
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        for (int k = 0; k < dim; k++) {
          for (int l = 0; l < dim; l++) {
            mat(qDitToDecimal(dim, {l,i}), qDitToDecimal(dim,{j,k})) = Aijkl(i,j,k,l);
          }
        }
      }
    }

    return mat;
  }

  std::vector<MatrixXd> computeSab(MatrixXd &FF, int maxdim){


    BDCSVD<MatrixXd> bdcsvd( FF,  ComputeThinU | ComputeThinV);
    //int info = bdcsvd.info();


    double maxsing = bdcsvd.singularValues().maxCoeff();

    //bdcsvd.setThreshold(pow(10, -10)/abs(maxsing));
    //int mi = bdcsvd.rank();
    int mi = bdcsvd.singularValues().size();
    if (mi > maxdim){
      mi = maxdim;
    }

    VectorXd sqrsingvec = bdcsvd.singularValues().cwiseSqrt();
    sqrsingvec = sqrsingvec(seq(0, mi-1)).eval();
    MatrixXd Sb = bdcsvd.matrixV()(all, seq(0, mi-1)).eval();
    MatrixXd Sa = bdcsvd.matrixU()(all,seq(0, mi-1)).eval();
    Sb = ( Sb*sqrsingvec.asDiagonal()).eval();

    Sa = (Sa*sqrsingvec.asDiagonal()).eval();

    std::vector<MatrixXd> stor = {Sa, Sb};


    return stor;
  }


  Tensor<double, 3> Convmat(MatrixXd &mat){
      int qnew = 0;
      int qold = 0;
      Tensor<double, 3> T;
      qnew = mat.cols();
      qold = int(sqrt(mat.rows()));
      T.resize(qold, qold, qnew);
      for (int i = 0; i < qnew; i++) {
        for (int k = 0; k < int(pow(qold, 2)); k++) {
          std::vector<int> vec = twodecToBaseQ(k, qold);
          int prvi = vec[0];
          int drugi = vec[1];
          T(prvi, drugi, i) = mat(k, i);
          }
        }
      return T;
      }

  Tensor<double, 3> StoT(MatrixXd &S){
    Tensor<double, 3> TS = Convmat(S);
    Eigen::array<IndexPair<int>, 1> fst = { IndexPair<int>(0, 1) };
    Tensor<double, 4> temp = TS.contract(TS, fst);
    Eigen::array<IndexPair<int>, 2> last = { IndexPair<int>(0, 0),IndexPair<int>(2, 1) };
    Tensor<double, 3> Temp2 = temp.contract(TS, last);
    auto arr = Temp2.dimensions();
    std::cout << "[" <<arr[0]<< ","<<arr[1]<<","<<arr[2]<<"]" <<'\n';

    return Temp2;
  }

  double Ancontr(Tensor<double,3> & Ta, Tensor<double,3> & Tb){
    Eigen::array<IndexPair<int>, 1> fst = { IndexPair<int>(1, 0)};
    Tensor<double, 4> T1 = Ta.contract(Tb, fst);

    Eigen::array<IndexPair<int>, 2> sec = { IndexPair<int>(1, 3),IndexPair<int>(2, 0) };
    Tensor<double, 4> T2 = T1.contract(T1, sec);
    Eigen::array<IndexPair<int>, 4> last = { IndexPair<int>(0, 2),IndexPair<int>(1, 1),
      IndexPair<int>(2, 3),  IndexPair<int>(3, 0)};
    Eigen::Tensor<double, 0> ContractedA = T2.contract(T1, last);
    double sum = ContractedA(0);

    return sum;
  }

  std::vector<Tensor<double,3>> RGstep(Tensor<double, 4> &Aijkl,int maxdim){

    MatrixXd mat = AtoM(Aijkl);

    std::vector<MatrixXd> Sab = computeSab(mat, maxdim);

    std::vector<Tensor<double, 3>> Tab = {StoT(Sab[0]), StoT(Sab[1])};
    //std::cout << Tab[0] << '\n';

    return Tab;

  }
  double calc_norm_factor(VectorXd &norms){
    double sum = 0;
    int N = norms.size();
    for (int i = 0; i < N; i++) {
      sum += pow(3,(N-i)) * log(norms(i));
    }
    return sum;
  }

  Tensor<double, 6> AtoQ(Tensor<double, 4> &A){
    Eigen::array<IndexPair<int>, 1> fst = { IndexPair<int>(0, 2) };
    Tensor<double ,6> temp = A.contract(A, fst);
    Eigen::array<IndexPair<int>, 2> last = { IndexPair<int>(1, 0),IndexPair<int>(3, 2)};
    Tensor<double, 6> temp2 = temp.contract(A, last);

    return temp2;

  }

  // ADD NORMALIZATION LATER!!:
  Vector2d Partitsum(int numq, double beta, double J, int steps, int maxdim){
    std::vector<Tensor<double, 3>> T0init = T0ab(numq,beta,J);
    Tensor<double, 3> Ta = T0init[0];
    Tensor<double, 3> Tb = T0init[1];



    Tensor<double, 4> A = Aijkl(Ta, Tb);
    VectorXd norms(steps);
    int ind = 0;
    while ((steps) > pow(10, -4)) {
      //Tensor<double, 6> TT = AtoQ(A);
      Eigen::Map<Eigen::VectorXd> vec(A.data(), A.size());

      double norm = vec.norm();
      norms(ind) = norm;
      A = A/norm;
      std::vector<Tensor<double,3>> Tab = RGstep(A, maxdim);
      Ta = Tab[0];
      Tb = Tab[1];

      A = Aijkl(Ta, Tb);

      steps -= 1;
      ind +=1;
    }

    double res = Ancontr(Ta, Tb);
    //double fact = calc_norm_factor(norms);
    double fact=0;
    Vector2d resall;
    resall << res, fact;
    return resall;
  }

  void lnZ_sweep_beta_tofile(int numq,double J, int steps, int maxdim, std::string imedat){
      VectorXd betasm = VectorXd::LinSpaced(39, 0.05, 1.0);
      VectorXd betalg = VectorXd::LinSpaced(120, 1.025, 4.0);
      VectorXd beta(betasm.size()+ betalg.size());
      beta <<betasm , betalg;
      VectorXd lnz(beta.size());
      for (size_t i = 0; i < beta.size(); i++) {
        std::cout << (1.0*i/(beta.size()-1))*100 << '\n';
        auto resall = Partitsum(numq, beta(i), J, steps, maxdim);
        lnz(i) = log(resall(0)) - resall(1);
      }

      std::ofstream outputFile(imedat, std::ios::out);
      outputFile << "beta, lnZ, db=" << beta(1)-beta(0) << "\n";
      for (size_t i = 0; i < beta.size(); i++) {
        outputFile << beta(i) << " " << lnz(i) << "\n";
      }

      outputFile.close();
    }



}

#endif
