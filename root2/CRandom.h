#ifndef random_
#define random_
#include <cstdlib>
#include <cmath>
#include "TRandom3.h"

/**
 * !\brief Random numbers for a number of distributions
 *
 * Random number generation using the root TRandom3 generator
 * this is an alternative to the C++ generator in the main GEMINI++ 
 * distribution
 */


class CRandom
{
 protected:
  CRandom();
  static  CRandom* fInstance; //!< instance member to make tis a singleton
  bool one; //!< used for Gaus
  float angle; //!< used for Gaus
  float x; //!< parameter
  float  pi; //!< 3.14159
  TRandom3 * randy;
 public:
  static CRandom* instance(); //!< instance member to make this a singleton

  ~CRandom();
  double Rndm();
  float Gaus(float mean,float sigma);
  float expDecayTime(float width);
  float BreitWigner(float mean ,float width);
};


#endif
