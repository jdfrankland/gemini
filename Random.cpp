#include "CRandom.h"



CRandom* CRandom::fInstance = 0;
/**
 * Makes an instance of CRandom
 */
CRandom* CRandom::instance() 
{
    if (fInstance == 0) {
        fInstance = new CRandom;
    }
    return fInstance;
}

/**
 * COnstructor for CRandom
 */
CRandom::CRandom()
{

  pi = acos(-1.);
  one = true;
  angle = 0.;
  x = 0.;

}

CRandom::~CRandom()
{

}

/**
 * Returns a random number with uniform distribution between 0 and 1
 */
double CRandom::Rndm()
{
  return  (float)rand()/(float)RAND_MAX;
}
//***********************
/**
 * Returns random number with Gaussian distribution
\param mean is the mean value of the Gaussian distribution
\param sigma is the standard deviation of the distribution
*/
float CRandom::Gaus(float mean, float sigma)
{
  float r;
  if (one)
  {
    x = sqrt(-2.*log(Rndm()+1.e-37));
     angle = 2.*pi*Rndm();
     one = 0;
     r = x*cos(angle);
    }
  else 
    {
      r = x*sin(angle);
      one = 1;
    }

  
  return sigma*r + mean;
}
//*************************
/**
 * returns a decay time in zs sampled from a exponential distribution
\param width is the total decay width in MeV
*/
float CRandom::expDecayTime(float width)
{
  // returns a time from an exponential decay distribution
  //consistent with the total decay width
  if(width>0.)
    return - 0.65824*log(Rndm()+1.0e-37)/width;
  else
    return 1.E30;
}
//****************************************
/**
 * Returns a random number with a BreitWigner distribution
\param mean is the mean value of the distribution
\param width is the Full Width Half max of the distribution
*/

float CRandom::BreitWigner(float mean, float width)
{
  float xx = (1.-2.*Rndm())*pi/2.;
  float yy = tan(xx)/2.;

  return yy*width + mean;
}
