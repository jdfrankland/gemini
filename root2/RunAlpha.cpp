//sample program to use GEMINI++
// in alpha induced reactions.
// The reaction cross section predicted by the 
// optical model with Avrigeanu's global parmeters
// is assumed to go completely into fusion with the spin distribution
// predicted by the optical model.
// The compound nuclei are then decayed with gemini++
// The xn cross sections are printed out.



#include "CRunAlpha.h"
#include "CNucleus.h"
#include "CAlphaOM.h"
#include "CMass.h"

CRunAlpha::CRunAlpha(int iZ, int iA, double Elab, int Nevents)
{
  CAlphaOM om((double)iZ,(double)iA,Elab);

  int iZCN = iZ + 2;
  int iACN = iA + 4;

  CNucleus CN(iZCN,iACN);
  CLevelDensity::setLittleA(12.);
  CN.setEvapMode(1);

  CN.printParameters();

  CMass * mass = CMass::instance();
  


  float Ecm = Elab*(float)iA/(float)(iACN);

  float Q = mass->getExpMass(iZ,iA) + 2.242 - mass->getExpMass(iZCN,iACN);
  cout << "Q = " << Q << " MeV" << endl;

  float Ex = Ecm + Q;
  cout << "Ex = " << Ex << " MeV" << endl;

  int xn[20]={0};
  int Nfission = 0;


  for (int i=0;i<Nevents;i++)
    {
      float fL = (float)om.getL(CN.ran->Rndm());
      CN.setCompoundNucleus(Ex,fL);
      CN.decay();
     if (CN.abortEvent)
       {
	 CN.reset();
         continue;
       }


     if (CN.isSymmetricFission()) Nfission++;
     else
       {
        int iStable = CN.getNumberOfProducts();

        CNucleus *productER = CN.getProducts(iStable-1);
	if (productER->iZ == iZCN)
	  {
            int x = iACN - productER->iA;
            if (x < 20) xn[x]++;
	  }
       }

   CN.reset();

  }
  cout << "tot x sec = " << om.getXsec() << " mb" << endl;

  float konst = om.getXsec()/(float)Nevents;

  cout << "sigFission = " << (float)Nfission*konst << " mb " << 
    sqrt((float)Nfission)*konst << endl;

  for (int i=0;i<10;i++)
    {
      cout << i << " " << (float)xn[i]*konst << " " <<
	sqrt((float)xn[i])*konst << endl;
    }

}
