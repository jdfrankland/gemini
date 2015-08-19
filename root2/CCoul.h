#ifndef coul_
#define coul_
#include <iostream>
#include <cmath>
#include <complex>
using namespace std;

class CCoul
{
 public:
  double logDerF(int ,double, double );
  complex<double> logDerH(int, double, double);
  double abSum (complex<double>);
  int init(int, double, double);
  double F;
  double dF;
  double G;
  double dG;
};
#endif
