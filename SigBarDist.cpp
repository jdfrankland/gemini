#include "CSigBarDist.h"

float CSigBarDist::width=1.;
float const CSigBarDist::width0=1.5;

//****************************************************************
  /**
   * constructor
   /param sName0 is the name of the files containing fitted coeff.
  */
CSigBarDist::CSigBarDist(string sName0, float Zp0, float Ap0)
{

  Zp = Zp0;
  Ap = Ap0;
  string sName = sName0;
  sigCharged[1] = new CSigCharged(sName,Zp,Ap);


  sName = sName0+"P";
  //see if second file is there

  string fullName=string(GINPUT)+"tl/"+sName+".inv";

  ifstream ifFile(fullName.c_str());

  if (ifFile.fail() || sName0 == "neutron" )
    {
      one = 1;
      return;
    }
  
  ifFile.close();
  ifFile.clear();
  one = 0;

  sigCharged[2] = new CSigCharged(sName,Zp,Ap); 
  sName = sName0+"M";
  sigCharged[0] = new CSigCharged(sName,Zp,Ap); 

}
//****************************************************
  /**
   * destructor
   */
CSigBarDist::~CSigBarDist()
{
  if (one) delete sigCharged[1];
  else for (int i=0;i<3;i++) delete sigCharged[i];
}

//*******************************************
/**
 * returns the quantity 
 * \f$S=\sum_{\ell=0}^{\infty} (2\ell+1)T_{\ell}(\varepsilon)\f$
 * which is related to the inverse cross section by
 * \f$S=\frac{\sigma_{inv}}{\pi\lambda^{2}}\f$
 \param fEk is the kinetic energy of the evaporated particle
 \param temp is temperature of daughter in MeV
 */
float CSigBarDist::getInverseXsec(float fEk, float temp)
{
  if (one || temp <= 0. || width == 0.) 
                          return sigCharged[1]->getInverseXsec(fEk);
  float deltaR = sqrt(temp)*width;
  float ee[3];
  for (int i=0;i<3;i++) 
    {
     ee[i] = sigCharged[i]->getInverseXsec(fEk);
     if (ee[i] == 0.) return 0.;
    }

   
  /*
  float c1 = (ee[2]-ee[0])/2./width0;
  float c2 = (ee[2]+ee[0]-2.*ee[1])/pow(width0,2)/2.;
  float s0 =  ee[1] + deltaR*c1 + c2*pow(deltaR,2);
  float s1 = ee[1];
  float s2 = ee[1] - deltaR*c1 + c2*pow(deltaR,2);
  return (s0+s1+s2)/3.;
  */
  float c2 = (ee[2]+ee[0]-2.*ee[1])/pow(width0,2)/2.;
  float  out =  ee[1] + 2./3.*c2*pow(deltaR,2);
  if (out < 0.) out = 0.;
  return out;
}
//**************************************************
/**
 * set the parameter controlling the width of the barrier distribution
\param width00 - radial shift is \f$ \Delta R= \sqrt T* width00 \f$
 */
void CSigBarDist::setBarWidth(float width00)
{
  width = width00;
}
//***************************************************
/**
 * returns the parameter controlling the width of the barrier dist
 */
float CSigBarDist::getBarWidth()
{
  return width;
}
//**************************************************************************
/**
 * prints out the width parameter
 */
void CSigBarDist::printParameters()
{
  cout << "tl barrier width parameter = " << width << endl;
}

//**************************************************************************
  /**
   *  prepares for a series of opertions for a given iZ
   /param iZ0 is proton number of daughter
   */
void CSigBarDist::prepare(float Z0, float A0)
{
  Z = Z0;
  A = A0;
  sigCharged[1]->prepare(Z,A);
  if (one) return;
  sigCharged[0]->prepare(Z,A);
  sigCharged[2]->prepare(Z,A);
}
//**********************************************************
  /**
   * returns the barrier in MeV
   */
float CSigBarDist::getBarrier()
{
  return sigCharged[1]->getBarrier();
}
