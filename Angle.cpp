#include "CAngle.h"

using namespace std;

// following line is needed in ROOT version
//ClassImp(CAngle)

float const CAngle::pi=acos(-1.);

/**
 * Constructor
 */
CAngle::CAngle(float theta0, float phi0)
{
  theta = theta0;
  phi = phi0;
}
//************************************************
CAngle CAngle::transform(CAngle angle1, CAngle angle2)
/**
 * This subroutine performs a rotational transformation.
 *
 * theta1,phi1 are the coordinates (spherical angles in
 * radians)of a unit vector in the original coordinate systems.
 * the z  axis is made to rotate in the phi=phi2 plane by an
 *angle theta2. The coordinates of the vector in the new
 *reference frame are theta3,phi3

 \param angle1 is the initial (theta,phi) angles 
 \param angle2 specifies the (theta,phi) rotation
 */

{

  //
  // To calculate the angle between two vectors
  // theta1,phi1 and theta2,phi2
  // use CALL transform(theta1,phi1,-theta2,phi2,theta3,phi3)
  // i.e. theta3 in the angle between them, note the negative sign
  // on theta2


  // rotate vector in x-y plane by -phi2
  // and find cartesian coordinates
  float xp = sin(angle1.theta)*cos(angle1.phi-angle2.phi);
  float yp = sin(angle1.theta)*sin(angle1.phi-angle2.phi);
  float zp = cos(angle1.theta);

  // rotate vector in z-x plane by theta2
  float zt = zp*cos(angle2.theta) - xp*sin(angle2.theta);
  float xt = xp*cos(angle2.theta) + zp*sin(angle2.theta);
  float yt = yp;

  // rotate vector in x-y plane back by +phi2
  // and find spherical coordinates
  float theta3, phi3;
  if (xt == 0.0 && yt == 0.0) phi3 = 0.;
  else
   {
     phi3 = atan2(yt,xt) + angle2.phi;
     if (phi3 >2.* pi) phi3 = phi3 - 2.*pi;
     if (phi3 < 0.0) phi3 = phi3 + 2.*pi;
   }

 if (zt > 1.0) theta3 = 0.;
 else if (zt < -1.0) theta3 = pi;
 else theta3 = acos(zt);

 CAngle angle3(theta3,phi3);
 return angle3;  
}
