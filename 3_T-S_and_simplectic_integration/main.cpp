#include "ts.hpp"
#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include <math.h>
#include <cmath>
# define PI           3.14159265358979323846
#include <fstream>
using std::ofstream;
#include "gnuplot-iostream.h"


int main(){

  double T = 20000;
  double dt = 0.01;
  Vector4d zac(1.0,0,0,0.5);


  osc Oscilator1(T, dt, 10, zac);
  osc Oscilator2(T, dt, 0.5, zac);
  osc Oscilator3(T, dt, 1, zac);
  osc Oscilator4(T, dt, 10, zac);

  iterate(Oscilator1, 0, "nezapisi");
  iterate(Oscilator2, 0, "ne");
  iterate(Oscilator3, 0, "ne");
  iterate(Oscilator4, 0, "ne");


  //std::cout << (1/Oscilator1.time(last))*(Oscilator1.storage(0,all).sum() - Oscilator1.storage(1, all).sum()) << '\n';
  //std::cout << abs(Oscilator1.Ear(last) - Oscilator1.E0) << '\n';
  //write_differences(Oscilator1, Oscilator2, Oscilator3, Oscilator4);
  return 0;
}
