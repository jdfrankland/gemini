#include "CRandom.h"
#include <cmath>

/**
 *!\brief angular distributions 
 * 
 * class to randomly selects a polar angle theta
 * of an emitted fragment
 *
 */


class CAngleDist
{
 protected:
  CRandom *ran; //!<random number generator
  float ** dist; //!<array containing sampled distributions

  static int const maxL; //!<maximum angular distribution for 
                         //!<which distributions are made


  static int const nAngle; //!< number of angular bins used 
                           //!<to sample distributions
  static float const pi; //!< the mathematical constant 
 public:

  CAngleDist();
  ~CAngleDist();


  float getTheta(int l);
};
