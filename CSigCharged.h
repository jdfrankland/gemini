#ifndef sigCharged_
#define sigCharged_

#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;

/**
 *!\brief inverse xsections
 *
 * This class calculates \f$ \sum_{0}^{\infty} (2\ell+1) T_{\ell}(\varepsilon)\f$ which is the inverse xsection divided by \f$ \pi/k^2 \f$
 * 
 */


class CSigCharged
{
 private:
  string sName; //!< name of input file of coefficients
  float rc0; //!< paramter to calculate radius for Coulomb barrier
  float rc1; //!< paramter to calculate radius for Coulomb barrier
  float rc2; //!< paramter to calculate radius for Coulomb barrier
  float omega0;  //!< paramter for omega
  float omega1;  //!< parameter for omega
  float omega2;  //!< parameter for omega
  float omega3;  //!< parameter for omega
  float rI0; //!< paramter for radius for rotational energy
  float rI1; //!< paramter for radius for rotational energy
  float rI2; //!< paramter for radius for rotational energy
  float aa0; //!< below barrier correction parameter
  float aa1; //!< below barrier correction parameter
  float a0; //!< above barrier correction parameter
  float a1; //!< above barrier correction parameter
  float Zp;  //!< proton number of evaporated particle
  float Ap;  //!< mass number of evaporated particle

  float barrier; //!< Coulomb barrier in MeV
  float InvInertia; //!< inverse of the moment of inertia associated with rotateion
  float omega; //!< amega parameter in MeV
  float a; //!< above barrier correction
  float aa; //!< below barrier correction
  float offset;  //!< offset

  bool neutron; //!<bool to signify neutron calculation
  float n0; //!< neutron parameter
  float n1; //!< neutron parameter
  float n2; //!< neutron parameter

 public:
  CSigCharged(string file, float Zp0, float Ap0); //!< constructor
  void prepare(float Z,float A); //!< prepares for calculations of inverse xsections
  float getInverseXsec(float energy); //!<calculates the inverse xsection 
  float getBarrier(); //!< calculates the barrier

 };
#endif
