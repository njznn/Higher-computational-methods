using namespace std;
#include <fstream>
using std::ofstream;
#include "dmrg.hpp"
#include "mps_algorithm.hpp"
#include <iomanip>
#include <chrono>

/*
MatrixXd entropybipart(std::string robni){
  VectorXi nspin(6);
  nspin << 4,6,8,10,12,14;


  SparseMatrix HeisH = makeSparseH<14>(robni);
  VectorXd base1 = basestatefromsparse(HeisH);
  auto dec = MPSdecomp<double>(base1,14, 2, "neodrezi", 0.0, 0.0);
  VectorXd maxen = entangeled_entropy<double>(dec[3]);
  MatrixXd entr = MatrixXd::Zero(maxen.size(),6);
  entr(all, 5) = maxen;

  for (size_t i = 0; i < 5; i++) {
    MatrixXd HeisenbergH = Heisenberg_H(nspin(i), robni);
    VectorXd base = basestate(HeisenbergH);
    auto decmop = MPSdecomp<double>(base,nspin(i), 2, "neodrezi", 0.0, 0.0);
    VectorXd res = entangeled_entropy<double>(decmop[3]);
    entr(seq(0, (res.size()-1)), i) = res;
  }

  return(entr);

}

MatrixXd entropybipartrand(){
  VectorXi nspin(6);
  nspin << 4,6,8,10,12,14;
  VectorXi N(6);
  N<< pow(2,4), pow(2,6), pow(2,8), pow(2,10), pow(2,12), pow(2,14);


  VectorXcd randomstatemax(N(5));
  randomstatemax.real() = VectorXd::Random(N(5));
  randomstatemax.imag() = VectorXd::Random(N(5));
  randomstatemax.normalize();
  auto decmop = MPSdecomp<complex<double>>(randomstatemax,nspin(5), 2, "neodrezi", 0.0, 0.0);
  VectorXd res = entangeled_entropy<complex<double>>(decmop[3]);
  MatrixXd entr = MatrixXd::Zero(res.size(),6);
  entr(all, 5) = res;


  for (size_t i = 0; i < 5; i++) {
    VectorXcd randomstate(N(i));
    randomstate.real() = VectorXd::Random(N(i));
    randomstate.imag() = VectorXd::Random(N(i));
    randomstate.normalize();
    auto decmop = MPSdecomp<complex<double>>(randomstate,nspin(i), 2, "neodrezi", 0.0, 0.0);
    VectorXd res = entangeled_entropy<complex<double>>(decmop[3]);
    entr(seq(0, (res.size()-1)), i) = res;
  }

  return(entr);

}


VectorXd checkMPSerr(std::string robni){
  VectorXd res(6);
  VectorXi nspin(6);
  nspin << 4,6,8,10,12,14;

  for (size_t i = 0; i < 5; i++) {
    MatrixXd HeisenbergH = Heisenberg_H(nspin(i), robni);
    VectorXd base = basestate(HeisenbergH);
    auto decmop = MPSdecomp<double>(base,nspin(i), 2, "neodrezi", 0.0, 0.0);
    if (nspin(i) ==4){
      VectorXd basedec = veckoeffrmps<double, 4>(decmop[0], 2);
      res(i) = (base-basedec).norm();
    }
    else if (nspin(i) ==6){
      VectorXd basedec = veckoeffrmps<double, 6>(decmop[0], 2);
      res(i) = (base-basedec).norm();
    }
    else if (nspin(i) ==8){
      VectorXd basedec = veckoeffrmps<double, 8>(decmop[0], 2);
      res(i) = (base-basedec).norm();
    }
    else if (nspin(i) ==10){
      VectorXd basedec = veckoeffrmps<double, 10>(decmop[0], 2);
      res(i) = (base-basedec).norm();
    }
    else if (nspin(i) ==12){
      VectorXd basedec = veckoeffrmps<double, 12>(decmop[0], 2);
      res(i) = (base-basedec).norm();
    }
  }
  SparseMatrix HeisH = makeSparseH<14>(robni);
  VectorXd base1 = basestatefromsparse(HeisH);
  auto dec = MPSdecomp<double>(base1,14, 2, "neodrezi", 0.0, 0.0);
  VectorXd basedec = veckoeffrmps<double,14>(dec[0], 2);
  res(5) = (base1-basedec).norm();

  return res;
}
*/
VectorXcd checkMPSerrrandom(VectorXcd & random, double tresh){
  VectorXcd res(6);
  VectorXi nspin(6);
  nspin << 4,6,8,10,12,16;
  VectorXi N(6);
  N<< pow(2,4), pow(2,6), pow(2,8), pow(2,10), pow(2,12), pow(2,14);



  for (size_t i = 4; i < 5; i++) {


    auto decmop = MPSdecomp<complex<double>>(random,nspin(i), 2, "odrezi", tresh);
    if (nspin(i) ==4){
      VectorXcd basedec = veckoeffrmps<complex<double>, 4>(decmop[0], 2);
      res(i) = (random-basedec).norm();
    }
    else if (nspin(i) ==6){
      VectorXcd basedec = veckoeffrmps<complex<double>, 6>(decmop[0], 2);
      res(i) = (random-basedec).norm();
    }
    else if (nspin(i) ==8){
      VectorXcd basedec = veckoeffrmps<complex<double>, 8>(decmop[0], 2);
      res(i) = (random-basedec).norm();
    }
    else if (nspin(i) ==10){
      VectorXcd basedec = veckoeffrmps<complex<double>, 10>(decmop[0], 2);
      res(i) = (random-basedec).norm();
    }
    else if (nspin(i) ==12){
      VectorXcd basedec = veckoeffrmps<complex<double>, 12>(decmop[0], 2);
      res(i) = (random-basedec).norm();
    }
    else if (nspin(i) ==14){
      VectorXcd basedec = veckoeffrmps<complex<double>, 14>(decmop[0], 2);
      res(i) = (random-basedec).norm();
  }


}
  return res;
}



int main(int argc, char const *argv[]) {






  //psi << 0.1,0.2,0.05,0.025,0.025,0.1,0.2,0.05,0.06,0.04,0.03,0.02,
  //0.03,0.02,0.02,0.03;






  std::cout << "break" << '\n';
  //std::cout << entropybipartrand() << '\n';
  /*
  fstream myfile;
  std::string ime_dat = "/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/8_DMRG_I/randvec.txt";
  myfile.open(ime_dat,fstream::out);
  */
  VectorXcd randomstate((int) pow(2,12));
  randomstate.real() = VectorXd::Random((int)pow(2,12));
  randomstate.imag() = VectorXd::Random((int)pow(2,12));
  /*
  for (size_t i = 0; i < randomstate.size(); i++) {
    myfile << randomstate.real()(i) << '\t' << randomstate.imag()(i) << '\n';
  }
  */

  std::cout << checkMPSerrrandom(randomstate, 0.001) << '\n';
  /*
  MatrixXd HeisenbergH4 = Heisenberg_H(4, "neperiodic");
  MatrixXd HeisenbergH6 = Heisenberg_H(6, "neperiodic");
  MatrixXd HeisenbergH8 = Heisenberg_H(8, "neperiodic");
  MatrixXd HeisenbergH10 = Heisenberg_H(10, "neperiodic");
  MatrixXd HeisenbergH12 = Heisenberg_H(12, "neperiodic");

  //SparseMatrix HeisH = makeSparseH<14>("periodic");

  VectorXd base4 = basestate(HeisenbergH4);

  VectorXd base6 = basestate(HeisenbergH6);
  VectorXd base8 = basestate(HeisenbergH8);
  VectorXd base10 = basestate(HeisenbergH10);
  VectorXd base12 = basestate(HeisenbergH12);

  //VectorXd base14 = basestatefromsparse(HeisH);
  */


  VectorXd entrab(6);
  VectorXi N(6);
  N<< pow(2,4), pow(2,6), pow(2,8), pow(2,10), pow(2,12), pow(2,14);
  std::cout << N << '\n';

  VectorXcd base4(N(0));
  base4.real() = VectorXd::Random(N(0));
  base4.imag() = VectorXd::Random(N(0));
  base4.normalize();
  VectorXcd base6(N(1));
  base6.real() = VectorXd::Random(N(1));
  base6.imag() = VectorXd::Random(N(1));
  base6.normalize();
  VectorXcd base8(N(2));
  base8.real() = VectorXd::Random(N(2));
  base8.imag() = VectorXd::Random(N(2));
  base8.normalize();
  VectorXcd base10(N(3));
  base10.real() = VectorXd::Random(N(3));
  base10.imag() = VectorXd::Random(N(3));
  base10.normalize();
  VectorXcd base12(N(4));
  base12.real() = VectorXd::Random(N(4));
  base12.imag() = VectorXd::Random(N(4));
  base12.normalize();


  VectorXcd base14(N(5));
  base14 = VectorXd::Random(N(5));
  base14.imag() = VectorXd::Random(N(5));
  base14.normalize();

/*
  entrab(0) = ABABentropysim<complex<double>,4>(base4, pow(10, -6));
  entrab(1) = ABABentropysim<complex<double>,6>(base6,pow(10, -6));
  entrab(2) = ABABentropysim<complex<double>,8>(base8,pow(10, -6));
  entrab(3) = ABABentropysim<complex<double>,10>(base10,pow(10, -6));
  entrab(4) = ABABentropysim<complex<double>,12>(base12,pow(10, -6));



  entrab(5) = ABABentropysim<complex<double>,14>(base14, pow(10, -6));
  std::cout << entrab << '\n';
  std::cout << "random" << '\n';
  */
  auto start = std::chrono::high_resolution_clock::now();


  //fstream myfile;
  //std::string ime_dat = "/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/8_DMRG_I/entropy_all_comp_pbc.txt";

  //myfile.open(ime_dat,fstream::out);
  //auto res = MPSdecomp(psi,4, 2, "neodrezi",2, 1.0);
  //std::cout << entangeled_entropy<4>(res[3]) << '\n';
  //std::cout << truncerr(res[2][0], res[3][0]) << '\n';
  //std::cout <<veckoeffrmps<4>(res[0],2) << '\n';
  /*
  VectorXd S = entropybipart("nonper");
  VectorXd SPBC = entropybipart("periodic");
  myfile << "nonper"<< "\t" << "pbc" << "\n";
  for (size_t i = 0; i < S.size(); i++) {
    myfile << S(i) << "\t" << SPBC(i) << "\n";
  }
  */
  /*
  VectorXd err = checkMPSerr("non");
  VectorXd errpbc = checkMPSerr("periodic");
  myfile << "nonper"<< "\t" << "pbc" << "\n";
  for (size_t i = 0; i < err.size(); i++) {
    myfile << err(i) << "\t" << errpbc(i) << "\n";
  }

  MatrixXd allbipart = entropybipart("periodic");
  myfile << "n=4"<< "\t"<<"::2" <<"\t"<< "n=14" << "\n";
  for (size_t i = 0; i < allbipart.rows(); i++) {
    for (size_t j = 0; j < allbipart.cols(); j++) {
      myfile<<allbipart(i,j)<< "\t";
    }
    myfile << "\n";
  }
  */


  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  cout << "Time taken by function: "
        << duration.count() << " microseconds" << endl;






  return 0;
}
