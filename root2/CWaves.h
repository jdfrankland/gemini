#include <iostream>

using namespace std;
class CWaves
{
 public:
  CWaves(double,double,int);
  ~CWaves();
  double  *F; // regular wavefunction
  double *dF; // derivative of regular
  double *G; // irregular wavefunction
  double *dG; // derivative of irregular
  double *sigma; // Coulomb phase shifts
 private:
  void coulombWave();
  void sphericalBessel();
  double rho;
  double gamma;
  int lMaximum;
};
//*******************************************************
