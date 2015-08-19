#ifndef random_
#define random_
#include <cstdlib>
#include <cmath>

/**
 * !\brief Random numbers for a number of distributions
 *
 * Random number generation using the C++ random number function
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
 public:
  static CRandom* instance(); //!< instance member to make this a singleton
  ~CRandom();
  double Rndm();
  float Gaus(float mean,float sigma);
  float expDecayTime(float width);
  float BreitWigner(float mean ,float width);
};


#endif
