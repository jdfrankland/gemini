#include "CMass.h"
//folowing is needed in ROOT version
//#include "Rtypes.h"

/**
 *!\brief fusion xsection, Bass model, excitation energy in fusion
 *
 * class to determine excitation energy and 
 *the critical angular momemtum in fusion.
 * The fusion xsection \f$ \sigma(\ell) = \pi\lambda^2 \sum \frac{2\ell +1}
{1+\exp\left(\frac{\ell-\ell_0}{\Delta_{\ell}}\right)}\f$
* where \f$\lambda\f$ is really lambdabar
* \f$ \Delta_{\ell}\f$ is the diffuseness
*/

class CFus
{
 protected:
  int iZp; //!< projectile proton number
  int iAp; //!<projectle mass number
  int iZt; //!<target proton number
  int iAt; //!<target mass number
  float fElab; //!< lab energy in MeV
  float R12; //!< sum of radii
  float U; //!< reduced mass
  float A; //!< const for Coulomb potential 
  float B; //!< const for centrifugal potential
  float C; //!< nuclear potential constant
  static float const D; //!< Bass potential parameter
  static float const E; //!< Bass potential parameter
  static float const G; //!< Bass potential parameter
  static float const H; //!< Bass potential parameter
  float E1; //!< critical energy 1 in Bass Model
  float E2; //!< critical energy 2 in Bass Model

  float MAX; //!< maximum L for fusion barrier

  float CL1; //!< angular momentum assocaited with E1
  float CL2; //!< angular momentum associated with E2
  float W[300]; //!< fusion barrier for each L

  float F(float R,float AL);
  float FF(float R,float AL);
  float  FFF(float R,float AL);
  float  FFFF(float R, float AL);
 public:
  float plb; //!< pi-lambdabar-squared in mb
  float dif; //!<diffuseness
  float Ecm; //!<reaction center of mass energy in MeV
  float vcm; //!<Compound nucleus velocity in cm/ns
  float vbeam; //!<beam velocity in cm/ns
  float Ex; //!< excitation energy
  int iZcn; //!< compound nucleus atomic number
  int iAcn;//!< compound nucleus mass number
  CFus(float plb0,float dif0);
  CFus(int iZprojectile, int iAprojectile, 
       int iZtarget, int iAtarget, float ELab, float dif0);
  void init(float plb0,float dif0);
  float getL0(float xsection);
  float getBassL(); 
  float getBassXsec(); 

  //following is needed in ROOT version
  //ClassDef(CFus,1); //Gemini CFus
};
