#ifndef QB_HPP
#define QB_HPP
#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include <math.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <bitset>
using std::ofstream;

static constexpr int N =12; // CHAIN SIZE! zaradi funkcije std::bitset je potrebno
// def N ze tukaj, lahko bi napisal svoje binarno-decimalne pretvorbe, potem bi
//v inicializaciji strukture dolocil velikost verige.

 // chain size
struct qb{ //qubit chain
  Vector<complex<double>, Dynamic> statekoef; //state
  Vector<complex<double>, Dynamic> statekoefcop; // need copy for autocorelation!
  // ni potrebno za ostale stvari!!
  //BETTER if i crate 2 chains with same initial random state!!

  complex<double> z; // can be -it for time or -\beta

  double dt; // cas korak
  double dbeta; // temp


  FILE * fout;

  qb( double dt_or_beta,  string str, VectorXd *ptr=NULL){
    srand((unsigned int) time(0));
    if (ptr !=NULL){
      statekoef = *ptr;
      statekoef.normalize();
      statekoefcop = statekoef;
    }
    else{
      statekoef = Matrix<complex<double>, int(pow(2,N)), 1> ::Random(pow(2, N),1);
      statekoef.normalize();
      statekoefcop = statekoef;

    }
    // specification what z is, time or beta:
    if (str=="cas"){
      z = (-1.0)*dt_or_beta* 1i;
      dt = dt_or_beta;
    }
    else if (str=="beta"){
      z = -dt_or_beta;
      dbeta = dt_or_beta;
    }


  }
  ~qb(){};
};

/// propagator U_i_i+1 on bits i_i+1 on state with coef ind
void U_i_oncoef(qb &ch, Vector<complex<double>, Dynamic> &copy,
  Vector<complex<double>, Dynamic> &copy2, unsigned int ind, int i , complex<double> a){
  auto bind = bitset<N>(ind);
  string twobits = (to_string(bind[(i+1)% N]) + to_string(bind[i]));


  if (twobits == "11" || twobits=="00"){
    ch.statekoef(ind) = exp(-1.0*ch.z * a)*exp(2.0*ch.z * a) * copy(ind);
    ch.statekoefcop(ind) = exp(-1.0*ch.z * a)*exp(2.0*ch.z * a) * copy2(ind);
  }
  else if (twobits == "01"){
    auto swapped = bind.flip(i).flip((i+1)%N);
    int swapind = stoi(swapped.to_string(), nullptr, 2);
    ch.statekoef(ind) = exp((-1.0)*ch.z*a)*(cosh(2.0*ch.z*a)*copy(ind)
    + sinh(2.0*ch.z*a)*copy(swapind));

    ch.statekoefcop(ind) = exp((-1.0)*ch.z*a)*(cosh(2.0*ch.z*a)*copy2(ind)
    + sinh(2.0*ch.z*a)*copy2(swapind));

  }

  else if(twobits == "10"){
    auto swapped = bind.flip(i).flip((i+1)%N);
    int swapind = stoi(swapped.to_string(), nullptr, 2);

    ch.statekoef(ind) = exp((-1.0)*ch.z*a)*(cosh(2.0*ch.z*a)*copy(ind)
    + sinh(2.0*ch.z*a)*copy(swapind));

    ch.statekoefcop(ind) = exp((-1.0)*ch.z*a)*(cosh(2.0*ch.z*a)*copy2(ind)
    + sinh(2.0*ch.z*a)*copy2(swapind));
  }

}


/// Propagator grabs i,i+1 bit and works on all states (vector of koef) at pair of bits i,i+1
void U_ij(qb & ch, int i, complex<double>a ){
  auto copy = ch.statekoef;
  auto copy2 = ch.statekoefcop;
  for (unsigned int ind = 0; ind < pow(2, N); ind++) {
    U_i_oncoef(ch,copy, copy2, ind, i, a);

  }
}


void U_even(qb & ch, complex<double>a){
  for (size_t i = 0; i < N/2; i++) {
    int i_ev = 2*i;
    U_ij(ch, i_ev, a);

  }
}

void U_odd(qb &ch, complex<double> a){
  for (size_t i = 0; i < N/2; i++) {
    int i_od = 2*i + 1;
    U_ij(ch, i_od, a);
  }
}

void Asym (qb &ch, double a){
	U_even (ch,a*1.0 + 0.0i);
  U_odd (ch,a*1.0 + 0.0i);

}

void S2_prep(qb &ch, double a){
  U_even(ch, a*0.5);
  U_odd(ch, a*1.0);
  U_even(ch, a*0.5);
}
void S2(qb &ch){
  U_even(ch, 0.5);
  U_odd(ch, 1.0);
  U_even(ch, 0.5);
}

// S4 is not working well, dont know why
void S4(qb & ch, double b){
  double f = pow(2, 1.0/3),
	       x0 = (-1)*f/(2.0 - f),
	       x1 = 1.0/(2 - f);


  S2_prep(ch, b*x1);
  S2_prep(ch, b*x0);
  S2_prep(ch, b*x1);
}

void S3_prep (qb & ch, int con, double a)
{
  complex<double> p1 = 0.125 * (1.0+ 1i/sqrt(3));
  complex<double> p2 = 2.0*p1;
  complex<double> p4 = conj(p2);
  complex<double> p5 = conj(p1);
  double p3 = 0.25;

	p1 *= a;
	p2 *= a;
	p3 *= a;
	p4 *= a;
	p5 *= a;

	if (con == 0)
	{
		U_even (ch, p1);
		U_odd (ch, p2);
		U_even (ch, p3);
		U_odd (ch, p4);
		U_even (ch, p5);
	}

	else if (con == 1)
	{
		U_even (ch, p5);
		U_odd (ch, p4);
		U_even (ch, p3);
		U_odd (ch, p2);
		U_even (ch, p1);
	}
}


void S3(qb &ch, double a){
  S3_prep(ch, 0, 2.0*a);
}

void S5(qb & ch, double a){
  S3_prep(ch, 0,1.0*a);
  S3_prep(ch, 1,1.0*a);
}

void sigmaz(qb &ch, int j, string w){  // sigma_z_j on j-th state, can work on copy or statekoef
  if (w=="normal"){
    for (unsigned int ind = 0; ind < pow(2, N); ind++) {
      auto bind = bitset<N>(ind);
      if (bind[j] == 1){
        ch.statekoef(ind) = (-1.0)* ch.statekoef(ind);
      }
    }
  }
  else if(w=="copy"){
    for (unsigned int ind = 0; ind < pow(2, N); ind++) {
      auto bind = bitset<N>(ind);
      if (bind[j] == 1){
        ch.statekoefcop(ind) = (-1.0)* ch.statekoefcop(ind);
      }
    }
  }
}

void propagate_time(qb &ch, string method, int final_t){
  int size = final_t/ch.dt;
  for (size_t i = 0; i < size; i++) {
    if (method=="S2"){
      S2(ch);
      std::cout << ch.statekoefcop.dot(ch.statekoefcop) << '\n';
    }
    else if (method=="S1"){
      Asym(ch,1.0);
      std::cout << ch.statekoef << '\n';
    }
    else if (method=="S3"){
      S3(ch, 1.0);
      std::cout << ch.statekoef.dot(ch.statekoef) << '\n';
    }
    else if(method=="S4"){
      S4(ch, 1.0);

    }
    else if (method=="S5"){
      S5(ch,1.0);
    }
  }
}

VectorXd sigmaz_autocol_one(qb &ch, int j, int final_t){
  int size = final_t/ch.dt;
  VectorXd Cj = VectorXd::Zero(size,1);
  Cj(0) = 1.0;
  for (size_t i = 1; i < size; i++) {
    sigmaz(ch, j, "normal");
    S2_prep(ch, 1.0);
    sigmaz(ch, j, "normal");
    Cj(i) = ch.statekoefcop.dot(ch.statekoef).real();
  }
  return Cj;
}

template<int objcount>
VectorXd Cz(qb (&obj)[objcount], int j , double final_t){
  qb ch = obj[0];
  int size = final_t/ch.dt;
  MatrixXd zbeta = MatrixXd::Zero(size,1);
  for (size_t i = 0; i < objcount; i++) {
    //zbeta(all, 0) += calc_ZorH_beta_one(obj[i], final_beta, "Asym",kol);
    zbeta(all, 0) += sigmaz_autocol_one(obj[i], j,final_t );
    std::cout << "delež narejenega:" << (double) (i+1)/objcount<< '\n';
    //zbeta(all, 2) += calc_ZorH_beta_one(obj[i], final_beta, "S3",kol);
    //zbeta(all, 3) += calc_ZorH_beta_one(obj[i], final_beta, "S4",kol);
    //zbeta(all, 4) += calc_ZorH_beta_one(obj[i], final_beta, "S5",kol);

  }

  zbeta = zbeta/objcount;
  return zbeta;

}

//spinski tok
void Jj_on_state(qb &ch, Vector<complex<double>, Dynamic> copy, int j){ //not working on copy!
  for (size_t ind = 0; ind < pow(2, N); ind++) {
    auto bind = bitset<N>(ind);
    string twobits = (to_string(bind[(j+1)% N]) + to_string(bind[j]));

    if (twobits == "11" || twobits=="00"){
      ch.statekoef(ind) = 0.0 + 0.0i;
      }
    else if (twobits == "01"){
      auto swapped = bind.flip(j).flip((j+1)%N);
      int swapind = stoi(swapped.to_string(), nullptr, 2);
      ch.statekoef(ind) = +2.0*1i*copy(swapind);

    }
    else if (twobits == "10"){
      auto swapped = bind.flip(j).flip((j+1)%N);
      int swapind = stoi(swapped.to_string(), nullptr, 2);
      ch.statekoef(ind) = -2.0*1i*copy(swapind);

    }
  }
}

  Vector<complex<double>, Dynamic> Jj_on_state_copy(qb &ch, int j){
    Vector<complex<double>, Dynamic> newstate =
     Matrix<complex<double>, int(pow(2,N)), 1> ::Zero(pow(2, N),1);

    for (size_t ind = 0; ind < pow(2, N); ind++) {
      auto bind = bitset<N>(ind);
      string twobits = (to_string(bind[(j+1)% N]) + to_string(bind[j]));

      if (twobits == "01"){
        auto swapped = bind.flip(j).flip((j+1)%N);
        int swapind = stoi(swapped.to_string(), nullptr, 2);
        newstate(ind) = 2.0*1i*ch.statekoef(swapind);

      }
      else if (twobits == "10"){
        auto swapped = bind.flip(j).flip((j+1)%N);
        int swapind = stoi(swapped.to_string(), nullptr, 2);
        newstate(ind) = -2.0*1i*ch.statekoef(swapind);
      }
  }

  return newstate;
}


VectorXd spin_current_one(qb &ch, int final_t){
  int size = final_t/ch.dt;
  VectorXd jcur = MatrixXd::Zero(size,1);
  auto copy = ch.statekoef;
  Jj_on_state(ch,copy, 0); // state chi!
  for (size_t t = 0; t < size; t++) {
    for (size_t n = 0; n < N; n++) {
      jcur[t] += ch.statekoefcop.dot(Jj_on_state_copy(ch, n)).real();
    }
    S2_prep(ch, 1.0);
  }
  return jcur;
}

template<int objcount>
VectorXd spin_current(qb (&obj)[objcount], double final_t){
  qb ch = obj[0];
  int size = final_t/ch.dt;
  VectorXd jcur = MatrixXd::Zero(size,1);
  for (size_t i = 0; i < objcount; i++) {
    jcur += spin_current_one(obj[i], final_t);
    std::cout << "delež narejenega:" << (double) (i+1)/objcount<< '\n';
  }
  return jcur/objcount;
}

// hamiltonian on coefficient (State vector) in :
Vector<complex<double>, Dynamic> H_on_state(qb &ch ){
  Vector<complex<double>, Dynamic> newvec = Matrix<complex<double>, Dynamic, 1> ::Zero(pow(2, N),1);
  for (size_t i = 0; i < N; i++) {
    for (unsigned int ind= 0; ind < pow(2, N); ind++) {
      auto bind = bitset<N>(ind);
      string twobits = (to_string(bind[(i+1)% N]) + to_string(bind[i]));
      if (twobits == "11" || twobits=="00"){
        newvec(ind) +=  ch.statekoef(ind);
      }
      else if (twobits == "01"){
        auto swapped = bind.flip(i).flip((i+1)%N);
        int swapind = stoi(swapped.to_string(), nullptr, 2);
        newvec(ind) += (-1.0)*ch.statekoef(ind)
        + 2.0* ch.statekoef(swapind);

      }

      else if(twobits == "10"){
        auto swapped = bind.flip(i).flip((i+1)%N);
        int swapind = stoi(swapped.to_string(), nullptr, 2);
        newvec(ind) += (-1.0)*ch.statekoef(ind)
        + 2.0* ch.statekoef(swapind);
      }
    }
  }
  return newvec;
}


/////// if a=0 return Z(beta) if a=1 return H'(Beta) (without 1/Z(beta))
VectorXd calc_ZorH_beta_one(qb ch, double final_beta, string method, int a){
  int size = final_beta/ch.dbeta;
  double db = ch.dbeta;
  VectorXd zbeta = VectorXd::Zero(size,1);
  VectorXd H = VectorXd::Zero(size,1);
  if (a==1){
    for (size_t i = 0; i < size; i++) {
        if (method=="Asym"){
          Asym(ch, 0.5);
          H(i) = ch.statekoef.dot(H_on_state(ch)).real();
        }
        else if(method=="S2"){
          S2_prep(ch, 0.5);
          H(i) = ch.statekoef.dot(H_on_state(ch)).real();
        }
        else if(method=="S4"){
          S4(ch, 0.5);
          H(i) = ch.statekoef.dot(H_on_state(ch)).real();

        }
        else if(method=="S3"){
          S3(ch, 0.5);
          H(i) = ch.statekoef.dot(H_on_state(ch)).real();

        }
        else if(method=="S5"){
          S5(ch, 0.5);
          H(i) = ch.statekoef.dot(H_on_state(ch)).real();
        }
      }
    return H;
  }
  else if(a==0){
    for (size_t i = 0; i < size; i++) {
      if (method=="Asym"){
        Asym(ch, 0.5);
        zbeta(i) =  ch.statekoef.dot(ch.statekoef).real();
      }
      else if(method=="S2"){
        S2_prep(ch, 0.5);
        zbeta(i) =  ch.statekoef.dot(ch.statekoef).real();
      }
      else if(method=="S4"){
        S4(ch, 0.5);
        zbeta(i) =  ch.statekoef.dot(ch.statekoef).real();

      }
      else if(method=="S3"){
        S3(ch, 0.5);
        zbeta(i) =  ch.statekoef.dot(ch.statekoef).real();

      }
      else if(method=="S5"){
        S5(ch, 0.5);
        zbeta(i) =  ch.statekoef.dot(ch.statekoef).real();
      }
    }
    return zbeta;
  }

  else {
    std::cout << "izberi pravi output" << '\n';
    exit(1);
    }
};



template<int objcount>
MatrixXd Zb(qb (&obj)[objcount], double final_beta, int kol){
  qb ch = obj[0];
  double db = ch.dbeta;
  int size = final_beta/ch.dbeta;
  MatrixXd zbeta = MatrixXd::Zero(size,5);
  for (size_t i = 0; i < objcount; i++) {
    zbeta(all, 0) += calc_ZorH_beta_one(obj[i], final_beta, "Asym",kol);
    zbeta(all, 1) += calc_ZorH_beta_one(obj[i], final_beta, "S2",kol);
    zbeta(all, 2) += calc_ZorH_beta_one(obj[i], final_beta, "S3",kol);
    zbeta(all, 3) += calc_ZorH_beta_one(obj[i], final_beta, "S4",kol);
    zbeta(all, 4) += calc_ZorH_beta_one(obj[i], final_beta, "S5",kol);
    std::cout << "delež narejenega:" << (double) (i+1)/objcount<< '\n';

  }

  zbeta = zbeta/objcount;
  return zbeta;
}


void write_mat_txt(qb &ch, MatrixXd mat, string ime_dat){
  int rows = mat.rows();
  int cols = mat.cols();

  fstream myfile;
  myfile.open(ime_dat,fstream::out);
  myfile << "dbeta = " << ch.dbeta<<" " <<"betamax = "<<rows*ch.dbeta << "\n";
  myfile << "Asym "<< "\t"<< "S2"<< "\t" << "S3"<< "\t" << "S4" << "\t" << "S5" << "\n";
  //myfile << "dt = " << ch.dt<<" " <<"tmax = "<<rows*ch.dt <<" "<< "S5" << "\n";
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      myfile << mat(i,j) << "\t";
    }
    myfile << "\n";
  }
  myfile.close();

}











#endif
