#include "CLevelDensity.h"

CLevelDensity* CLevelDensity::fInstance = 0;

bool  CLevelDensity::normal = true;

float const CLevelDensity::pi=acos(-1.);
float CLevelDensity::k0=7.3;
float CLevelDensity::kInfinity=12.;
float CLevelDensity::aKappa = 0.00517;
float CLevelDensity::cKappa = .0345;
float CLevelDensity::af_an = 1.036;
float CLevelDensity::aimf_an = 1.02;
float CLevelDensity::eFade = 18.52;
float CLevelDensity::jFade = 50.;
float CLevelDensity::Ucrit0 = 9.;
float CLevelDensity::Jcrit = 14;
//float  CLevelDensity::Ucrit0 = 0.;
//float  CLevelDensity::Jcrit = -1.;

/**
 * Constructor
 */
CLevelDensity::CLevelDensity()
{
  //constructor read in in level density parameter

  
  string fileName("tbl/gemini.inp");
  string fullName;
  if (getenv("GINPUT") == NULL) fullName = fileName;
  else
    {
      string dir(getenv("GINPUT"));
     fullName = dir+fileName;
    }
  ifstream ifFile (fullName.c_str());

  if (ifFile.fail())
    {
      cout << "Gemini++: unable to open gemini.inp" << endl;
      return;
    }

  string line;
  string text1;
  string text2;
  string text3;
  string text4;
  getline(ifFile,line);
  ifFile >> text1 >> k0 >> text2 >> aKappa >> text3 >> cKappa >> text4
         >> kInfinity;

  ifFile >> text1 >> eFade >> text2 >> jFade >> text3;
  ifFile >> text1 >> Ucrit0 >> text2 >> Jcrit >> text3;

  ifFile.close();
  ifFile.clear();

  
}

CLevelDensity* CLevelDensity::instance() // mod-TU
{
    if (fInstance == 0) {
        fInstance = new CLevelDensity;
    }
    return fInstance;
}

//*********************************************************************
  /**
   * Returns the backshifted excitation energy energy to be used in the 
   * Fermi gas formula for the level Density
     \param fU0 is the thermal excitation energy in MeV
     \param fPairing is the pairing correction in MeV
     \param fShell is the shell correction in MeV
     \param fJ is the angular momentum in units of hbar
  */
float CLevelDensity::getU(float fU0, float fPairing, float fShell, float fJ)
{




  fU = fU0;
  if (fU <= 0.)
    {
      fU = 0.;
     return fU;
    }
  //simple fade out of pairing
  float shiftP;
  float Ucrit = 0.;
  if (fJ < Jcrit) Ucrit = Ucrit0*pow(1.-fJ/Jcrit,2);
  if (fU > Ucrit) shiftP = fPairing;
  else shiftP = fPairing*(1.-pow(1.-fU/Ucrit,2));

  fU += shiftP;
  if (fU <= 0.)
    {
      fU = 0.;
     return fU;
    }   
  //fU += fShell*(1.0-exp(-fU/eFade-fJ/jFade));
  float shiftS = fShell*tanh(fU/eFade+fJ/jFade);




  fU += shiftS;

  if (fU < 0.) fU = 0.;

  J = fJ;
  return fU;
}



//*********************************************
/**
 * Returns the level-density parameter in units of MeV-1, the getU function
 * must already have been called to use this version
\param iA is the mass number
\param fJ is the angular momentum in units of hbar
\param iFission is a short indicating we are detail with a saddle-point shape
 */

float CLevelDensity::getLittleA(int iA, short iFission/*=0*/)
{
 //calculates the level density parameter
  //iA is nucleus mass number
  //fU is thermal excitation energy
  //fPairing is the pairing energy
  //fShell is the shell correction to the mass

  float fA = (float)iA;
  float kappa = 0.;
  float daden_dU;
  if ((normal && fU/fA < 3.) || aKappa == 0.)
     {
       if (k0 == kInfinity)
	 {
           aden = fA/k0;
	   daden_dU = 0.;
	 }
       else
	 {
	   if (aKappa > 0.) kappa = aKappa*exp(cKappa*fA);
	   //kappa = 1.5+.1143*J;
          float expTerm = exp(-kappa*fU/fA/(kInfinity-k0));
          aden = fA/(kInfinity - (kInfinity-k0)*expTerm);
          daden_dU = -pow(aden/fA,2)*kappa*expTerm;
	 }
       switch(iFission)
	 {
	 case 1:
           aden *= af_an;
           daden_dU *= af_an;
	   break;
         case 2:
           aden *= aimf_an;
           daden_dU *= aimf_an;
	   break;
	 }
  }
  else 
  { 
     float r;
     if (iFission == 1)r = 1.0696;
     if (iFission == 2)r = 1.05;
     else r = 1.0;
     if (aKappa > 0.) kappa = aKappa*exp(cKappa*fA);
     float fofr = ( kInfinity - ( kInfinity - k0 ) * r ) / (k0*r);
     float expTerm = exp(-fofr*kappa*fU/fA/(kInfinity-k0));
     aden = fA/(kInfinity - (kInfinity-k0)*r*expTerm);
     daden_dU = -pow(aden/fA,2)*kappa*expTerm*fofr;
  }

  entropy = 2.*sqrt(aden*fU);
  temp = sqrt(fU/aden)/(1.+fU/aden*daden_dU);

  return aden;
}

//*********************************************
/**
 * Returns the level-density parameter in units of MeV-1
\param iA is the mass number
\param fU0 is the thermal excitation energy in MeV
\param fPairing is the pairing correction in MeV
\param fShell is the shell correction in MeV
\param fJ is the angular momentum in units of hbar
\param iFission is a short indicating we are detail with a saddle-point shape
 */

float CLevelDensity::getLittleA(int iA, float fU0, float fPairing/*=0.*/,
	float fShell/*=0.*/, float fJ/*=0.*/, short iFission/*=0*/)
{
 //calculates the level density parameter
  //iA is nucleus mass number
  //fU is thermal excitation energy
  //fPairing is the pairing energy
  //fShell is the shell correction to the mass


  //if (getU(fU0, fPairing,fShell,fJ) <= 0.) return 0.;
  getU(fU0, fPairing,fShell,fJ);
  return getLittleA(iA,iFission);
}

//************************************************
/**
 * Returns the natural logrithm of the spin dependent level density in MeV-1
 * if fJ < 0, then littleA is calculated taking into account the spin reduction
 * of paring, but subsequentally the level density is calculated for 
 * zero spin. This is done when using Weisskopf formalism for evaporation.
 \param iA is the mass number
 \param fU0 is the thermal excitation energy in MeV
 \param fPairing is the pairing energy in MeV
 \param fShell is the shell correction in MeV
 \param fJ is the spin in hbar
 \param fMinertia is the moment of inertia in nucleon mass* fm^2
\param iFission indicates its for the saddle-point configuration
 */

//spin dependent Fermi_gas level density
float CLevelDensity::getLogLevelDensitySpherical
  (int iA, float fU0, float fPairing,
   float fShell, float fJ, float fMinertia, short iFission/*=0*/)
{
  //calculates the level density
  //iA is nucleus mass number
  //fU is thermal excitation energy
  //fPairing is the pairing energy
  //fShell is the shell correction to the mass

  if (getLittleA(iA,fU0,fPairing,fShell,fabs(fJ),iFission) == 0.) return 0.;
  if (fU <=0.) return 0.;
  if (fJ < 0.) fJ = 0.;
  float sigma =  fMinertia*temp/40.848;
  float preExp = (2.*fJ+1.)/(1.+pow(fU,(float)(1.25))*pow(sigma,(float)1.5))
    /pow(aden,(float)0.25)/24./sqrt(2.);
  return entropy + log(preExp);
}
//******************************************
/**
 * Returns the temperature in MeV
 * getLittleA must be called first
 */
float CLevelDensity::getTemp()
{
  return temp;
}
//*************************************
  /**
   * Returns the entropy.
   * getLittleA must be called first
   */
float CLevelDensity::getEntropy()
{
  return entropy;
}
//*************************************
/**
* Returns the level-density parameter in MeV-1
 */
float CLevelDensity::getAden()
{
  return aden;
}
//************************************************
  /**
   * Returns the spin independent Fermi-gas level density
   \param iA is the mass number
   \param fU0 is the thermal excitation energy in MeV
   \param fPairing is the pairing correction in MeV
   \param fShell is the shell correction in MeV
   */
float CLevelDensity::getLogLevelDensitySpherical
  (int iA, float fU0, float fPairing, float fShell)
{
  //calculates the level density
  //iA is nucleus mass number
  //fU is thermal excitation energy
  //fPairing is the pairing energy
  //fShell is the shell correction to the mass

  if (getLittleA(iA,fU0,fPairing,fShell)== 0.) return 0.;
  if (fU <=0.) return 0.;

  float preExp = sqrt(pi)/12./pow(aden,(float)0.25)/
    (1.+pow(fU+temp,(float)1.25));
  return entropy + log(preExp);
}
//***************************************************
/**
 * Sets the constants for  level-density parameter
 * where \f$a=\frac{A}{k_{\infty} - \left(k_{\infty} -k_{0} \right) \exp\left( \frac{\kappa}{k_{\infty}-k_{0}}\frac{U}{A}\right)}\f$
 * where \f$ \kappa = a_{\kappa} \exp\left(c_{\kappa} A\right) \f$.
 * Note setLittleA(8.) is equivalent to \f$a=A/8\f$
\param k00 is \f$k_{0}\f$
\param aKappa0 is \f$a_{\kappa}\f$
\param cKappa0 is \f$c_{kappa}\f$
\param kInfinity0 is \f$k_{\infty}\f$
*/
void CLevelDensity::setLittleA(float k00, float aKappa0/*=0.*/,
  float cKappa0 /*=0.*/, float kInfinity0 /*=12.*/)
{
  k0 = k00;
  aKappa = aKappa0;
  cKappa = cKappa0;
  kInfinity = kInfinity0;
  if (aKappa == 0.) kInfinity = k0;
}
//****************************************************
 /**
 * returns the inverse level-density parameter at zero excitation energy
 */
float CLevelDensity::getK0()
{
  return k0;
}
//************************************************************
 /**
 * returns the inverse level-density parameter at infinite excitation
 * energy
 */
float CLevelDensity::getKInfinity()
{
  return kInfinity;
}
//****************************************************
 /**
 * returns one of the coefficients used to calculate kappa
 * \f$ \kappa = a_{\kappa} \exp\left(c_{\kappa} A\right) \f$
 */
float CLevelDensity::getAKappa()
{
  return aKappa;
}
//****************************************************
 /**
 * returns the other coefficients used to calculate kappa
 * \f$ \kappa = a_{\kappa} \exp\left(c_{\kappa} A\right) \f$

 */
float CLevelDensity::getCKappa()
{
  return cKappa;
}
//*****************************************************
/**
 * returns the ration of saddle-point to equilibrium level-density parameter
 * for symmetyric fission
 */
float CLevelDensity::getAfAn()
{
  return af_an;
}
//*****************************************************
/**
 * returns the ration of saddle-point to equilibrium level-density parameter
 * for asymmetric fission
 */
float CLevelDensity::getAimfAn()
{
  return aimf_an;
}
//*****************************************************
 /**
  * set the ration of level-density paramters at the saddle-point to
  * equilibrium shape for symmetric fission
  \param af_an0 is the level-density parameter ratio for saddle-point to equilibirum deformation
  */
void CLevelDensity::setAfAn(float af_an0)
{
  af_an = af_an0;
}
//*****************************************************
 /**
  * set the ration of level-density paramters at the saddle-point to
  * equilibrium shape for asymmetric fisison , ie. imf emission
  \param aimf_an0 is the level-density parameter ratio for saddle-point to equilibirum deformation
  */
void CLevelDensity::setAimfAn(float aimf_an0)
{
  aimf_an = aimf_an0;
}
//***************************************************
  /**
   * writes out the values of the parameters
   */
void CLevelDensity::printParameters()
{
  cout << "k0 = " << k0 << " kInfinity= " << kInfinity << endl;
  cout << "aKappa= " << aKappa << " cKappa= " << cKappa << endl;
  cout << "af/an= " << af_an << " for symmetric fission" <<endl;
  cout << "aimf/an= " << aimf_an << " for asymmetric fission" <<endl;
  cout << "eFade= " << eFade << " jFade= " << jFade << endl;
  cout << "Ucrit0= " << Ucrit0 << " Jcrit= " << Jcrit << endl;
}
//***************************************************
  /**
   * returns the level density for a scission configuration.
   * This is used for evaporation from saddle to scission and
   * for determining the fission mass asymmetry
\param iA is mass number of fissioning system
\param fU is thermal excitation of fission system
\param adenInv in the inverse level-density paramter in MeV, is \f$a=A/adenInv\f$
  */
float CLevelDensity::getLogLevelDensityScission(int iA, float fU,
						float adenInv/*=8.*/ )
{
  aden = (float)iA/adenInv;
  temp = sqrt(fU/aden);
  entropy = 2.*sqrt(aden*fU);

 float preExp = sqrt(pi)/12./pow(aden,(float)0.25)/
    (1.+pow(fU+temp,(float)1.25));
  return entropy + log(preExp);
}
//*****************************************************
 /**
  * set the critical thermal excitation where pairing vanishes
  \param Ucrit0 is critical thermal excitaion energy in MeV
  \param Jcrit  is critical angular momentum where pairing vanishes 
  */
void CLevelDensity::setUcrit(float Ucrit00, float Jcrit0)
{
  Ucrit0 = Ucrit00;
  Jcrit = Jcrit0;

}
