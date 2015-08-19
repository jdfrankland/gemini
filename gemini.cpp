#include "CNucleus.h"

/**
 * These external functions allow GEMINI to be run form fortran
 */


int const m=100;

extern "C"
{

  extern struct
  {
    int z[m];
    int a[m]; 
    float T[m];
    float P[m];
    float theta[m];
    float phi[m];
  } gem_;


  /**
   *  initialization of GEMINI. Does not need to be run if you are
   * using the default parameters
   */
  void geminiinit_()
  {
    CNucleus::setEvapMode(1);  //force Hauser-Feshbach
    //CLevelDensity::setAfAn(1.01);
  }
  //*********************************************************************
  /**
   * returns the decay width of a compound nucleus in MeV
  \param iZ is the proton number of the compound nucleus
  \param iA is the mass number of the compound nucleus
  \param fEx is the excitation energy of the compound nucleus
  \param fJ is the spin of the compound nucleus
   */
  void decaywidth_(int* iZ, int* iA, float *fEx, float*fJ, float* width,
     float* logleveldensity)
  {
    CNucleus CN(*iZ,*iA,*fEx,*fJ);
    *width = CN.getDecayWidth();
    *logleveldensity = CN.getLogLevelDensity();
  }
  //********************************************************************
  /**
   * Calculates the decay of the nucleus
\param iZ is the proton number of the compound nucleus
\param iA is the mass number of the compound nucleus
\param fEx is the excitation energy of the compound nucleus
\param fJ is the spin of the compound nucleus
  */


  int gemini_(int* iZ, int* iA, float *fEx, float *fJ, float* thetaJ, 
	      float* phiJ,float *vx, float *vy, float *vz )
  {

    CNucleus CN(*iZ,*iA,*fEx,*fJ);
    CN.setVelocityCartesian(*vx,*vy,*vz);
    CAngle angle(*thetaJ,*phiJ);
    CN.setSpinAxisDegrees(angle);
    CN.decay();

    int n;
    //check to see if gemini was successfully able to decay the fragment
    if (CN.abortEvent)
      {
	n = 0;
      }
    else
      {
    n= CN.getNumberOfProducts();
    if (n > m)
      {
	cout << " increase common block size" << endl;
        abort();
      }
    CNucleus* product;
   
    for (int i=0;i<n;i++)
      {
        product = CN.getProducts(i);
        gem_.z[i] = product->iZ;
        gem_.a[i] = product->iA;
        gem_.T[i]= product->getKE();
        gem_.P[i] = product->getMomentum();
        angle = product->getAngleDegrees();
        gem_.theta[i] = angle.theta;
        gem_.phi[i] = angle.phi;
      }
      }
    CN.reset();
    return n;
  }


}


