#ifndef _levelDensity
#define _levelDensity

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;

/**
 *!\brief returns level density, entropy, temperature
 *
 * This class deals with level densities and related paramters.
 * the level-density parameter is \f$a=\frac{A}{k_{\infty} - \left(k_{\infty} -k_{0} \right) \exp\left( \frac{\kappa}{k_{\infty}-k_{0}}\frac{U}{A}\right)}\f$
 * where \f$ \kappa = a_{\kappa} \exp\left(c_{\kappa} A\right) \f$
 */

class CLevelDensity
{
 private:
  CLevelDensity();
  static CLevelDensity *fInstance; //!< instance member to make this class a singleton
  static float k0; //!< inverse level-density parameter at U=0
  static float kInfinity; //!< inverse level-density parameter at U=infinity
  static float aKappa; //!< mass dependence of kappa
  static float cKappa; //!< mass dependence of kappa
  static float af_an; //!< ratio of sym-fission saddle  to equilibrium level-density para
  static float aimf_an; //!< ratio of asy-fission(imf) saddle to equilibrium level-density para
  static bool normal;
  static float Ucrit0; //!< Ucrit for J= 0, vanishing of pairing
  static float Jcrit; //!< spin for the vanishing of pairing
   float aden; //!< little-density parameter
  float entropy; //!< entropy
  float temp;  //!< temperature
  static float eFade; //!< fade out of shell effects with excitation energy
  static float jFade; //!< fade out of shell effects with spin
  float fU;   //!< thermal excitation energy in MeV
  static float const pi; //!< the mathematical constant \f$\pi\f$

 public:
  static CLevelDensity *instance(); //!< instance member to make this class a singleton
  float getLittleA(int iA,float fU0,float fPairing=0.,float fShell=0.,
                 float fJ=0.,short iFission=0);
  float getLittleA(int iA, short iFission);
  float getU(float fU0,float fPairing, float fShell, float fJ);
  float getAden();
  float getLogLevelDensitySpherical(int iA,float fU0,float fPairing,
	       float fShell,float fJ,float fMinertia,
               short iFission=0);
  float getLogLevelDensitySpherical(int iA,float fU0,float fPairing,
         float fShell);
  float getTemp();
  float getEntropy();
  static void setLittleA(float k00,float aKappa0=0., float cKappa0=0.,
          float kInfinity=12.);
  static void setAfAn(float af_an0);
  static void setAimfAn(float aimf_an0);
  static void setUcrit(float Ucrit00, float Jcrit);
  static float getAKappa();
  static float getCKappa();
  static float getK0();
  static float getKInfinity();
  static float getAfAn();
  static float getAimfAn();
  static void printParameters();
  float getLogLevelDensityScission(int iA, float U, float adenInv=8.);
  float J;
};
#endif
