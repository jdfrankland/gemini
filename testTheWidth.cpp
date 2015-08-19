#include <iostream>
#include "CNucleus.h"
#include <cmath>


using namespace std;


int main()
{
  int iZ = 12; 
  int iA = 24;
  float fJ = 0.;  //compound nucleus spin
  cout << "spin = " << fJ << endl;
  CNucleus CN(iZ,iA);
  CN.setEvapMode(1);  //set Hauser feshbach calculation

  for (int i=0;i<10;i++)
    {
      float fEx = 14. + (float)i; //compound nucleus excitation energy
      CN.excite(fEx,fJ);

      float width = CN.getDecayWidth(); //decay width in MeV 
      float logRho = CN.getLogLevelDensity();
      float rho = exp(logRho); //level density in MeV-1

      cout << "Ex = " << fEx << " width = " << width << " MeV, rho = " <<
	rho << " MeV-1, width*rho = " << width*rho << endl;
    }

}

