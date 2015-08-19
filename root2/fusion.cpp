#include "CRun.h"
#include "CFus.h"


int main()
{
  int numTot = 8000; //spectra
  float const d0=4.;  //diffuseness of CN spin distribution
  string title0("fusion"); //name of output root file without extension


  float Elab = 9.17*12.;  // bombarding energy of projectile in the lab frame
  
  //define projectile
  int iZp = 6;
  int iAp = 12;

  //define target
  int iZt = 28;
  int iAt = 58;

  //use the Bass model to find critical angular momentum

  CFus fus(iZp,iAp,iZt,iAt,Elab,d0);
  float l0 = fus.getBassL();


  cout << "Elab= " << Elab << " l0= " << l0
       << " plb = " << fus.plb << " mb " << endl;
  cout << "iZ= " << fus.iZcn << " iAcn= " << fus.iAcn 
   << " Ex= " << fus.Ex << endl;

  //run gemini
  CRun run(fus.iZcn,fus.iAcn,fus.Ex,l0,
  d0,(int)(l0+4.*d0),fus.plb,numTot,title0);

}
