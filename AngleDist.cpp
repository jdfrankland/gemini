#include "CAngleDist.h"

int const CAngleDist::maxL = 20;
int const CAngleDist::nAngle = 90;
float const CAngleDist::pi= acos(-1.);

//**********************************************************
/**
 * Constructor
 */
CAngleDist::CAngleDist()
{
  ran = CRandom::instance();
  dist = new float* [maxL];
  for (int i=0;i<maxL;i++)
    {
      dist[i] = new float [nAngle];

      for (int j=0;j<nAngle;j++)
	{
          float theta = (float)(j+1)/180.*pi;
          dist[i][j] = pow(sin(theta),2*i+1);
          if (j > 0) dist[i][j] += dist[i][j-1];
	}

      for (int j=0;j<nAngle;j++)
	{
          dist[i][j] /= dist[i][nAngle-1];
	}

    }


}
//***********************************************************
/**
 * Destructor
 */
CAngleDist::~CAngleDist()
{
 for (int i=0;i<maxL;i++)
   {
     delete [] dist[i];
   }
 delete [] dist;
}
//*************************************************************
//!returns a random polar angle associated with a specified angular momentum
/**
 * The angle is chosen from the
 * distribution \f$P_l^l[cos(\theta)]^2 sin(\theta)\f$.
 * This angle is expressed in radians.
 * If the angular momentum is above maxL, an angle of 
 * 90 degrees is returned
 \param l is the orbital angular-momentum quantum number
 */
float CAngleDist::getTheta(int l)
{
  if (l >=maxL) return pi/2.;

  float xran = ran->Rndm();
  int i=0;
  for (;;)
    {
      if (xran < dist[l][i]) break;
      i++;
      if (i == nAngle) break;
    }
  float theta = ((float)(i) + ran->Rndm())*pi/180.;
  if (ran->Rndm() > 0.5) theta = pi - theta;

  return theta;
}
