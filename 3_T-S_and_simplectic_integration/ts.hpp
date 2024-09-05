#ifndef TS_HPP
#define TS_HPP
#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include <math.h>
#include <cmath>
# define PI           3.14159265358979323846
#include <fstream>
using std::ofstream;
#include <fstream>


struct osc {

  double Tm, dt, L, N;
  Vector<complex<double>, 4> state;

  Matrix<double, 4, Dynamic> storage;
  int n; //indeks strenutne iteracije
  double E, E0;
  VectorXd Ear;
  VectorXd time;
  FILE *fout;

  osc(double T_, double dt_, double L_, Vector4d in_state):
  Tm(T_), dt(dt_), state(in_state), L(L_){
    state = in_state;
    E0 = 0.5*pow(real(state(2)),2) +  0.5*pow(real(state(3)),2) +
    0.5 *pow(real(state(0)),2) + 0.5*pow(real(state(1)),2) +
    L*pow(real(state(0)), 2) *pow(real(state(1)),2);
    N = Tm/dt;
    time = VectorXd::Zero(N,1);
    storage = MatrixXd::Zero(4, N);
    storage(all, 0) = in_state;
    Ear = MatrixXd::Zero(N, 1);
    Ear(0) = E0;
    n = 0;

  }
   ~osc(){};
};

template<typename B>
void U_T(osc & oscl, B c){

  oscl.state(0)  =  oscl.state(0) + oscl.state(2)*oscl.dt * c ;
  oscl.state(1)  =  oscl.state(1) + oscl.state(3)*oscl.dt * c ;

}

template< typename B>
void U_V(osc & oscl, B c){

  oscl.state(2) = oscl.state(2) - c*oscl.dt*oscl.state(0)*(1.0+2.0*oscl.L*pow(oscl.state(1),2));
  oscl.state(3) = oscl.state(3) - c*oscl.dt*oscl.state(1)*(1.0+2.0*oscl.L*pow(oscl.state(0),2));
}


void hamilt(osc & oscl){

  double T = 0.5*pow(real(oscl.state(2)),2) +  0.5*pow(real(oscl.state(3)),2);
  double V = 0.5 *pow(real(oscl.state(0)),2) + 0.5*pow(real(oscl.state(1)),2) +
  oscl.L * pow(real(oscl.state(0)),2)*pow(real(oscl.state(1)),2);
  oscl.E = T+V;
}


void S2(double c,  osc & oscl){
  U_T< double>(oscl,  c*0.5);
  U_V<double>(oscl, c*1.0);
  U_T<double>(oscl, c*0.5);
}


void S4(osc & oscl){
  double x0 = (-1.0)* pow (2, 1.0/3)/(2 - pow (2, 1.0/3)),
	       x1 = 1.0 / (2 - pow (2, 1.0/3));

  S2 (x1, oscl);
  S2 (x0, oscl);
  S2 (x1, oscl);
}


void S3(osc &oscl){
  complex<double> p1 = 0.25 * (1.0+ 1i/sqrt(3));
  complex<double> p2 = 2.0*p1;
  double p3 = 0.5;

  U_T<complex<double>>(oscl, conj(p1));
  U_V<complex<double>>(oscl, conj(p2));
  U_T<complex<double>>(oscl, p3);
  U_V<complex<double>>(oscl, p2);
  U_T<complex<double>>(oscl, p1);


}

void S3_bar(osc &oscl){
  complex<double> p1 = 0.25 * (1.0+ 1i/sqrt(3));
  complex<double> p2 = 2.0*p1;
  double p3 = 0.5;

  U_T<complex<double>>(oscl, p1);
  U_V<complex<double>>(oscl, p2);
  U_T<complex<double>>(oscl, p3);
  U_V<complex<double>>(oscl, conj(p2));
  U_T<complex<double>>(oscl, conj(p1));

}

void S3_all(osc & oscl){
  S3_bar(oscl);
  S3(oscl);

}

void shrani(osc &oscl){
  oscl.Ear(oscl.n) = oscl.E;
  oscl.storage(all, oscl.n) = oscl.state.real();
  oscl.time(oscl.n) = oscl.n*oscl.dt;
}

void dump(osc &oscl){
  fprintf(oscl.fout, "% 15lf \n", oscl.E);

}

void iterate(osc & oscl, int s, string out){
  if (out == "zapisi") {
    oscl.fout = fopen("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/3_T-S_and_simplectic_integration/rez.txt", "w");

    while(oscl.n < (oscl.N-1)){
      switch(s){

          case 1:
              S2(1.0,oscl);
              break;
          case 2:
              S4(oscl);
              break;
          case 3:
              S3_all(oscl);
              break;
          default:
              S4(oscl);
              break;
      }
      oscl.n++;
      hamilt(oscl);
      shrani(oscl);
      dump(oscl);


    delete (oscl.fout);
  }

}
  else{
    while(oscl.n < (oscl.N-1)){
      switch(s){

          case 1:
              S2(1.0,oscl);
              break;
          case 2:
              S4(oscl);
              break;
          case 3:
              S3_all(oscl);
              break;
          default:
              S4(oscl);
              break;
      }
      oscl.n++;
      hamilt(oscl);
      shrani(oscl);

  }
}
}

void write_differences(osc &oscl1,osc & oscl2,osc &oscl3, osc &oscl4){
  ofstream file;
  file.open("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/3_T-S_and_simplectic_integration/razlike.txt", std::ios_base::app);
  file << "dt "<<"l=0" << "l=0.5" << "l=1" << "l=10" << "\n";
  for (size_t i = 0; i < oscl1.N; i++) {
    file <<oscl1.time[i] << ' '<<abs(oscl1.Ear(i)- oscl1.E0) << ' ' << abs(oscl2.Ear(i)-oscl1.E0) << ' ' << abs(oscl3.Ear(i)-oscl3.E0) << ' ' <<
    abs(oscl4.Ear(i)-oscl4.E0) << "\n" ;
  }
  file.close();
}


#endif
