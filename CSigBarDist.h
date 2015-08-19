#include "CSigCharged.h"
#include <cstdlib>
/**
 *!\brief inverse cross section with barrier distributions
 * 
 * calculates transmission coefficientinverse cross sections 
 * using a simplistic barrier 
 * distribution logic. The final result is the average of three  
 * inverse cross sections.
 *  \f$ T_{l}(Ek) = \frac{T_{l}^{R_{0}-\Delta R}(Ek) + T_{l}^{R_{0}}(Ek) + T_{l}^{R_{0}+\Delta R}(Ek)}{3} \f$ 
 * where \f$ T_{l}^{R_{0}}\f$ is the standard coeff derived using the IWBC and 
 * global optical-model potential. The other two are for when the radius of 
 * the nuclear potential is shifted by \f$ \pm \Delta R \f$. 
 * where \f$  \Delta R = width*\sqrt{temperature} \f$.
 */


class CSigBarDist
{
 private:
  CSigCharged* sigCharged[3]; //!< arrays for standard radii and +- width0
  bool one; //!< if true, no distribution, just standard radius is used
  static float width; //!< width paramter determines shifted radii
  static float const width0; //!< results readin from file for this shift
  float Z; //!< calls to getInverseXsec refer to this residual proton number
  float A; //!< calls to getInverseXsec refer to this residual mass number
  float Zp; //!<proton number of evaporated particle
  float Ap; //!<mass number of evaporated particle
 public:
  CSigBarDist(string, float Zp0, float Ap0 );
  ~CSigBarDist();
  float getInverseXsec(float fEk, float temp);
  static void setBarWidth(float width00);
  static float getBarWidth();
  static void printParameters();
  void prepare(float Z, float A);
  float getBarrier();
};
