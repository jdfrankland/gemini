#include "CTlArray.h"

/**
 *!\brief transmission coefficients with barrier distibution
 * 
 * calculates transmission coefficients using a simplistic barrier 
 * distribution logic. The final result is the average of three transmission 
 * coefficients. 
 *  \f$ T_{l}(Ek) = \frac{T_{l}^{R_{0}-\Delta R}(Ek) + T_{l}^{R_{0}}(Ek) + T_{l}^{R_{0}+\Delta R}(Ek)}{3} \f$ 
 * where \f$ T_{l}^{R_{0}}\f$ is the standard coeff derived using the IWBC and 
 * global optical-model potential. The other two are for when the radius of 
 * the nuclear potential is shifted by \f$ \pm \Delta R \f$. 
 * where \f$  \Delta R = width*\sqrt{temperature} \f$.
 */


class CTlBarDist
{
 private:
  CTlArray* tlArray[3]; //!< arrays for standard radii and +- width0
  bool one; //!< if true, no distribution, just standard radius is used
  static float width; //!< width paramter determines shifted radii
  static float const width0; //!< results readin from file for this shift
  int iZ; //!< calls to getTl refer to this residual proton number

 public:
  CTlBarDist(string);
  ~CTlBarDist();
  float getTl(int iL,float fEk, float temp);
  float getTlLow(int iL,float fEk, float temp);
  float getTlHigh(int iL,float fEk, float temp);
  float getInverseXsec(float fEk, float temp);
  static void setBarWidth(float width00);
  static float getBarWidth();
  static void printParameters();
  void prepare(int iZ0);
};
