#include "CNucleus.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

/**
 *!\brief example of fusion with root histograms
 *
 * class that I use to simulate statistical decay in fusion 
 * reaction. It produces a number of histograms, such a mass, charge
 * and evaporation spectra - The output is stored in root file
 */ 


class CRun
{

 public:
CRun(int iZ, int iA, float fEx, float l0, float d0, int lmax, float plb,
     int numTot,string title0,float vcm=0.,float thetaDetMin=0.,
     float thetaDetMax = 360.);
};
