#include "CNucleus.h"
// this is an example of using GEMINI CNucleus class to give the
//statistical decay of a compound nucleus

int main()
{
  int iZCN = 82; // proton number of compound nucleus
  int iACN = 198; // mass number of compound nucleus
  CNucleus CN(iZCN,iACN); //constructor
  float fEx = 67.; //excitation energy of compound nucleus
  float fJ = 40; // spin of compound nucleus

  CN.setCompoundNucleus(fEx,fJ); //specify the excitation energy and spin

  CN.setVelocityCartesian(); // set initial CN velocity to zero
  CAngle spin(CNucleus::pi/2,(float)0.); 
  CN.setSpinAxis(spin); //set the direction of the CN spin vector

  for (int i=0;i<30;i++)
    {
     cout << "event = " << i << endl;
     //CN.setWeightIMF();// turn on enhanced IMF emission
     CN.decay(); //decay the compound nucleus
  
     if (CN.abortEvent)
       {
	 CN.reset();
         continue;
       }

     // print of number of stable 
     cout << "number of products= " <<CN.getNumberOfProducts() << endl; 
                                               // fragments produced in decay

     cout << "number of gamma rays= " <<CN.getnGammaRays() << endl;
                                               // gamma rays produced in decay
     cout << "sum of energy  of gamma rays= " <<CN.getSumGammaEnergy() << endl;

     CNucleus * products = CN.getProducts(0); //set pointer to firt
                                                    //stable product

     for(int i =0; i<CN.getnGammaRays(); i++)
       cout << i+1 << " gamma ray energy = " << CN.getGammaRayEnergy(i)<< endl;

     for(;;)
       {
	 CNucleus * parent;
         parent = products->getParent();
         if (parent == NULL) cout << "stable fragment= " 
          << products->getName() << endl;
	 else 
         cout << "stable fragment= " 
	      << products->getName() << 
	  " parent = " << parent->getName() << endl;
         products = CN.getProducts();  // go to the next product
         if (products == NULL) break;
       }

   CN.reset();
    }


}
