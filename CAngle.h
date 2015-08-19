// -*- mode: c++ -*- 
//
#ifndef angle_
#define angle_

#include <cmath>

//******ROOT*********
  //#include "Rtypes.h"

/**
 *!\brief polar angles
 *
 * Class to deal with polar angles
 */

class CAngle
{
protected:

 public:
  static float const pi; //!< 3.14159
  float theta; //!< polar angle in radians
  float phi; //!< azimuth angle in radians
  CAngle(float,float);
  CAngle(){};
  static CAngle transform(CAngle angle1,CAngle angle2);

  //*****ROOT************
    //ClassDef(CAngle,1) //Gemini CAngle
};


#endif
