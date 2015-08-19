#include "CNucleus.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

/**
 *!\brief example of fusion with root histograms
 *
 * class that I use to simulate statistical decay in fusion 
 * reactions where there is a thick target. 
 * It produces a number of histograms, such a mass, charge
 * and evaporation spectra - The output is stored in root file
 */ 


class CRunThick
{

 public:
  CRunThick(int iZ, int iA, float fEx_min,float fEx_max, float l0_min, 
       float l0_Max,float d0, int lmax, float plb,int nBins,
     int numTot,string title0,float vcm=0.,float thetaDetMin=0.,
     float thetaDetMax = 360.);
};
