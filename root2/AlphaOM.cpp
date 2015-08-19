#include "CAlphaOM.h"

using namespace std;


CAlphaOM::CAlphaOM(double Z, double A, double Elab)
{




  double Ecm = Elab*A/(A+4.);

  alpha.findPara(Z,A,Elab);
  scatter.init(alpha.mu,alpha.V,0.,alpha.R,alpha.a,alpha.zz,
	       alpha.Rc,alpha.Wvol,alpha.Wsur,alpha.Ri,alpha.ai);


   int l= 0;
   xsec = 0.;
   for (;;)
     {
      double tl = scatter.getTl_OM(Ecm,l);
      sl[l] = (double)(2*l+1)*tl*scatter.plb;
      xsec += sl[l];
      if (tl < .001) break;
      l++;
     }
   lmax  = l;


   for (int l=0;l<=lmax;l++)
     {
       prob[l] = sl[l]/xsec;
       if (l > 0) prob[l] += prob[l-1];
     }

}
//****************************************************************
  /**
   * returns the reaction cross section in mb for a given ell wave
   \param l ell wave number
   */
double CAlphaOM::getSigL(int l)
{
  if (l > lmax) return 0.;
  return sl[l];
}
//*******************************************************************
double CAlphaOM::getXsec()
{
  return xsec;
}
//***************************************************************
int CAlphaOM::getLmax()
{
  return lmax;
}
//************************************************************
int CAlphaOM::getL(double random)
{
  int l= 0;
  for (;;)
    {
      if (random < prob[l]) break;
      if (l == lmax ) break;
      l++;
    }
  return l;
}
