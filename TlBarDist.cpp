#include "CTlBarDist.h"

float CTlBarDist::width=1.;
float const CTlBarDist::width0=1.5;

//****************************************************************
  /**
   * constructor
   /param sName0 is the name of the files containing fitted coeff.
  */
CTlBarDist::CTlBarDist(string sName0)
{
  string sName = sName0;
  tlArray[1] = new CTlArray(sName);


  sName = sName0+"P";
  //see if second file is there 

  string fullName;
  if (getenv("GINPUT") == NULL) fullName = "tl/"+sName+".tl";
  else
    {
      string dir(getenv("GINPUT"));
      fullName = dir+"tl/"+sName+".tl";
    }
  ifstream ifFile(fullName.c_str());


  //do not include fluctuations for protons and neutrons
  //The IWBC barriers do not extrapolate very well for 
  // high spins when the barrier is increases and we often get
  //some artifacts, so its best not to include these.
  if (ifFile.fail() || sName0 == "neutron" || sName0 == "proton" )
    {
      one = 1;
      return;
    }
  
  ifFile.close();
  ifFile.clear();
  one = 0;

  tlArray[2] = new CTlArray(sName); 
  sName = sName0+"M";
  tlArray[0] = new CTlArray(sName); 

}
//****************************************************
  /**
   * destructor
   */
CTlBarDist::~CTlBarDist()
{
  if (one) delete tlArray[1];
  else for (int i=0;i<3;i++) delete tlArray[i];
}
//******************************************************
  /**
   * returns the transmission coeff, including barrier distibution
/param iL is orbital angular momentum of evaporated particle
/param fEk is the kinetic energy in MeV of the evaporated particle
/param temp is the temperature in MeV of daughter
  */ 
float CTlBarDist::getTl(int iL, float fEk, float temp)
{
  if (one || temp <= 0. || width == 0.) return tlArray[1]->getTl(iL,fEk);
  if (tlArray[0]->iZMin >= iZ) return tlArray[1]->getTl(iL,fEk);
  if (tlArray[2]->iZMin >= iZ) return tlArray[1]->getTl(iL,fEk);

  float deltaR = sqrt(temp)*width;
  float ee[3];
  for (int i=0;i<3;i++) ee[i] = tlArray[i]->getTermInExp(iL,fEk);


  // for proton emission light nuclei at high angular momentum, 
  // the parameterized 
  // transmission coefficients are extrapolations. For tlArray[0]
  // we sometime find this extrapolation is bad, the tl value is 
  // larger than for tlArray[1]. tlArray[2] extrapolates better
  // which is good as it defined the subbarrier behavoir
  // the following if block is desigend to correct for bad extrapolations
  //of tlArray[0]

  
  if (ee[0] < ee[1]) 
    {
      double fE = fEk-1.;
      double tlNew;
      double tlOld = ee[2];
      for(;;)
        {
          if (fE <= 0.) break;
          tlNew = tlArray[2]->getTermInExp(iL,fE);
          if (tlNew >= ee[1])break;
          fE--;
          tlOld = tlNew; 
        }
      if (fE <= 0.)
        {
         ee[0] = 1000.;
         if (ee[1] > ee[0]) ee[0] = ee[1];
        }
     else
       {
         fE -= (ee[1]-tlOld)/(tlNew-tlOld) -1.;
         ee[0] = tlArray[1]->getTermInExp(iL,fE);
       }

    }
  

  float c1 = (ee[2]-ee[0])/2./width0;
  float c2 = (ee[2]+ee[0]-2.*ee[1])/pow(width0,2)/2.;


  float Tl = 1./(1.+exp(ee[1]));
  //float eee = ee[1] + deltaR*(c1 + c2*deltaR);
  float eee = ee[1] + deltaR*c1 + c2*pow(deltaR,2);
  Tl += 1./(1.+exp(eee));
  //eee = ee[1] - deltaR*(c1 - c2*deltaR);
  eee = ee[1] - deltaR*c1 + c2*pow(deltaR,2);
  Tl += 1./(1.+exp(eee));
  Tl /=3.;
  return Tl;

}

//******************************************************
  /**
   * The transmission coeff is determine from the average of three
   * transmission coeff. This routine returns the one of these three
   * transmission coeff. with the lowest barrier
/param iL is orbital angular momentum of evaporated particle
/param fEk is the kinetic energy in MeV of the evaporated particle
/param temp is the temperature in MeV of daughter
  */ 
float CTlBarDist::getTlLow(int iL, float fEk, float temp)
{
  if (one || temp <= 0. || width == 0.) return tlArray[1]->getTl(iL,fEk);
  if (tlArray[0]->iZMin >= iZ) return tlArray[1]->getTl(iL,fEk);
  if (tlArray[2]->iZMin >= iZ) return tlArray[1]->getTl(iL,fEk);

  float deltaR = sqrt(temp)*width;
  float ee[3];
  for (int i=0;i<3;i++) ee[i] = tlArray[i]->getTermInExp(iL,fEk);


  // for proton emission light nuclei at high angular momentum, 
  // the parameterized 
  // transmission coefficients are extrapolations. For tlArray[0]
  // we sometime find this extrapolation is bad, the tl value is 
  // larger than for tlArray[1]. tlArray[2] extrapolates better
  // which is good as it defined the subbarrier behavoir
  // the following if block is desigend to correct for bad extrapolations
  //of tlArray[0]

  
  if (ee[0] < ee[1]) 
    {
      double fE = fEk-1.;
      double tlNew;
      double tlOld = ee[2];
      for(;;)
        {
          if (fE <= 0.) break;
          tlNew = tlArray[2]->getTermInExp(iL,fE);
          if (tlNew >= ee[1])break;
          fE--;
          tlOld = tlNew; 
        }
      if (fE <= 0.)
        {
         ee[0] = 1000.;
         if (ee[1] > ee[0]) ee[0] = ee[1];
        }
     else
       {
         fE -= (ee[1]-tlOld)/(tlNew-tlOld) -1.;
         ee[0] = tlArray[1]->getTermInExp(iL,fE);
       }

    }
  

  float c1 = (ee[2]-ee[0])/2./width0;
  float c2 = (ee[2]+ee[0]-2.*ee[1])/pow(width0,2)/2.;


  float eee = ee[1] + deltaR*c1 + c2*pow(deltaR,2);
  float Tl = 1./(1.+exp(eee));
  return Tl;

}
//******************************************************
  /**
   * The transmission coeff is determine from the average of three
   * transmission coeff. This routine returns the one of these three
   * transmission coeff. with the highest barrier

/param iL is orbital angular momentum of evaporated particle
/param fEk is the kinetic energy in MeV of the evaporated particle
/param temp is the temperature in MeV of daughter
  */ 
float CTlBarDist::getTlHigh(int iL, float fEk, float temp)
{
  if (one || temp <= 0. || width == 0.) return tlArray[1]->getTl(iL,fEk);
  if (tlArray[0]->iZMin >= iZ) return tlArray[1]->getTl(iL,fEk);
  if (tlArray[2]->iZMin >= iZ) return tlArray[1]->getTl(iL,fEk);

  float deltaR = sqrt(temp)*width;
  float ee[3];
  for (int i=0;i<3;i++) ee[i] = tlArray[i]->getTermInExp(iL,fEk);


  // for proton emission light nuclei at high angular momentum, 
  // the parameterized 
  // transmission coefficients are extrapolations. For tlArray[0]
  // we sometime find this extrapolation is bad, the tl value is 
  // larger than for tlArray[1]. tlArray[2] extrapolates better
  // which is good as it defined the subbarrier behavoir
  // the following if block is desigend to correct for bad extrapolations
  //of tlArray[0]

  
  if (ee[0] < ee[1]) 
    {
      double fE = fEk-1.;
      double tlNew;
      double tlOld = ee[2];
      for(;;)
        {
          if (fE <= 0.) break;
          tlNew = tlArray[2]->getTermInExp(iL,fE);
          if (tlNew >= ee[1])break;
          fE--;
          tlOld = tlNew; 
        }
      if (fE <= 0.)
        {
         ee[0] = 1000.;
         if (ee[1] > ee[0]) ee[0] = ee[1];
        }
     else
       {
         fE -= (ee[1]-tlOld)/(tlNew-tlOld) -1.;
         ee[0] = tlArray[1]->getTermInExp(iL,fE);
       }

    }
  

  float c1 = (ee[2]-ee[0])/2./width0;
  float c2 = (ee[2]+ee[0]-2.*ee[1])/pow(width0,2)/2.;


  float eee = ee[1] - deltaR*c1 + c2*pow(deltaR,2);
  float Tl = 1./(1.+exp(eee));
  return Tl;

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
float CTlBarDist::getInverseXsec(float fEk, float temp)
{
  float tot = 0.;
  float xmax = 0.;
  int iL = 0;
  for(;;)
    {
      float x = (float)(2*iL+1)*getTl(iL,fEk,temp);
      xmax = max(x,xmax);
      tot += x;
      if (x < xmax*.01) break;
      iL++;
      if (iL > 20) break;
    }
  return tot;
}
//**************************************************
/**
 * set the parameter controlling the width of the barrier distribution
\param width00 - radial shift is \f$ \Delta R= \sqrt T* width00 \f$
 */
void CTlBarDist::setBarWidth(float width00)
{
  width = width00;
}
//***************************************************
/**
 * returns the parameter controlling the width of the barrier dist
 */
float CTlBarDist::getBarWidth()
{
  return width;
}
//**************************************************************************
/**
 * prints out the width parameter
 */
void CTlBarDist::printParameters()
{
  cout << "tl barrier width parameter = " << width << endl;
}

//**************************************************************************
  /**
   *  prepares for a series of opertions for a given iZ
   /param iZ0 is proton number of daughter
   */
void CTlBarDist::prepare(int iZ0)
{
  iZ = iZ0;
  tlArray[1]->prepare(iZ);
  if (one) return;
  tlArray[0]->prepare(iZ);
  tlArray[2]->prepare(iZ);
}
