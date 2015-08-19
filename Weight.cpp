#include "CWeight.h"

float const fact0=0.3;
//****************************************************************************
/**
 * Determines the degree of weighting for enhance IMF emission
    \param gammaLight decay width for light-particle evaporation (MeV)
    \param gammaImf decay width for IMF emission (MeV)
    \param gammaFission Fission decay width (MeV)
    \param gammaGamma Gamma-ray decay width

*/
void CWeight::findFactor(float gammaLight, float gammaImf, float gammaFission,
float gammaGamma)
{
  if (gammaImf <= 0.) 
    {
      fact = 1.;
      iWeight = 0;
      return;
    }
  if (gammaLight+gammaFission+gammaGamma <=0.)
    {
      fact = 1.;
      iWeight = 0;
      return;
    }

  float gammaTot = gammaLight+gammaImf+gammaFission+gammaGamma;
  fact = fact0*(gammaTot-gammaImf)/(1.-fact0)/gammaImf;
  if (fact < 1.)fact = 1.; 
}
//********************************************************* 
  /**
   * Chooses the decay channel given the decay widths.
   * If IMF weighting is truned on, then IMF emission is enhanced and 
   * a running weight is calculated.
    \param gammaLight decay width for light-particle evaporation (MeV)
    \param gammaImf decay width for IMF emission (MeV)
    \param gammaFission Fission decay width (MeV)
    \param gammaGamma Gamma-ray decay width
    \param xran Random number
  */
    
int  CWeight::chooseChannel(float gammaLight, float gammaImf, 
			    float gammaFission, float gammaGamma, float xran)
{

  if (iWeight == 1) findFactor(gammaLight,gammaImf,gammaFission,gammaGamma);

  float fact1 = 1.;
  if (iWeight) fact1 = fact;  //weighting turned on
  

  float Gamma_real = gammaLight + gammaImf + gammaFission + gammaGamma;
  float Gamma_weight = gammaLight + gammaImf*fact1 + gammaFission + gammaGamma;
  
  float probLight = gammaLight/Gamma_weight;
  float probImf = gammaImf*fact1/Gamma_weight + probLight;
  float probFission = gammaFission/Gamma_weight + probImf;

  //no need of weighting if only imf emission
  if (gammaLight + gammaFission + gammaGamma <=0.) Gamma_weight = Gamma_real;


  int iChan;
  if (xran < probLight) 
    {
     iChan = 0;
     runningWeight *= Gamma_weight/Gamma_real;
     if (iWeight > 0)iWeight++;
    }
  else if (xran < probImf) 
    {
      iChan = 1;
      runningWeight *= Gamma_weight/fact1/Gamma_real;
      iWeight = 0; //turn weighting off for rest of event
    }
  else if (xran < probFission)
    {
      iChan =2;
      runningWeight *= Gamma_weight/Gamma_real;
      iWeight = 0; //turn weighting off for the rest of the event
     }
  else 
    {
      iChan =3;  //gamma ray emission
      runningWeight *= Gamma_weight/Gamma_real;
      iWeight = 0; //turn weighting off for the rest of the event
    }
 
  return iChan;

}
//************************************************************************
  /**
   * turns on IMF weighting, i.e. emhanced IMF emissions
   */
void CWeight::setWeightIMF()
{
  iWeight = 1;
}
//**********************************************************************
  /**
   * When IMF weighting is turn on, this routine returns the weighting factor
   * which should be used to increment all histograms. If weighting is turned
   * off, then unity is returned
   */
float CWeight::getWeightFactor()
{
  return runningWeight;
}
