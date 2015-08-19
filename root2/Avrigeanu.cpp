#include "CAvrigeanu.h"
#include <cmath>


//**************************************************
  /**
   * gives the global optical model potential for 
   * alpha particle scattering described by Avrigneau
   * PRC 49 (1994) 2136
   \param Z target Z
   \param A target A
   \param E lab energy of alpha particle
   */

void CAvrigeanu::findPara(double Z, double A, double E)
{
  double A3 = pow(A,1./3.);

  R = 1.245*A3;
  a = 0.817 - 0.0085*A3;
  V = 101.1 + 6.051*Z/A3 - 0.248*E;


  Ri = 1.57*A3;
  ai = 0.692-0.02*A3;
  if (E < 73.) Wvol = 12.64 - 1.706*A3 + 0.2*E;
  else Wvol = 26.82 - 1.706*A3 + .006*E;

  Rc = 1.17*A3;


  mu = A*4./(A+4.);
  zz = 2.*Z;

}
