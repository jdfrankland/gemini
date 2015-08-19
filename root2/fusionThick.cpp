#include "CRunThick.h"
#include "CFus.h"


int main()
{
  int numTot = 80000; //spectra
  float const d0=4.;  //diffuseness of CN spin distribution
  string title0("fusionThick"); //name of output root file without extension


  float Elab_max = 9.17*12.;  // bombarding energy at front or target
  float Elab_min = 8.74*12.;  // bombarding energy at back of target
  //the difference between these two is the 
  //energy loss of the projectile in the target
  
  //define projectile
  int iZp = 6;
  int iAp = 12;

  //define target
  int iZt = 28;
  int iAt = 58;

  //use the Bass model to find critical angular momentum

  //find critical spin for min bombarding energy
  CFus fus_min(iZp,iAp,iZt,iAt,Elab_min,d0);
  float l0_min = fus_min.getBassL();

  //find critical spin for maximum bombaring energy
  CFus fus_max(iZp,iAp,iZt,iAt,Elab_max,d0);
  float l0_max = fus_max.getBassL();


  cout << "Elab_min= " << Elab_min << " l0_min= " << l0_min
       << " plb_min = " << fus_min.plb << " mb " << endl;
  cout << "iZ= " << fus_min.iZcn << " iAcn= " << fus_min.iAcn 
   << " Ex_min= " << fus_min.Ex << endl;

  cout << "Elab_max= " << Elab_max << " l0_max= " << l0_max
       << " plb_max = " << fus_max.plb << " mb " << endl;
  cout << "iZ= " << fus_max.iZcn << " iAcn= " << fus_max.iAcn 
   << " Ex_max= " << fus_max.Ex << endl;



  //run gemini
  CRunThick run(fus_min.iZcn,fus_min.iAcn,fus_min.Ex,fus_max.Ex,l0_min,l0_max,
  d0,(int)(l0_max+4.*d0),
	   (fus_min.plb+fus_max.plb)/2.,10,numTot,title0);

}
