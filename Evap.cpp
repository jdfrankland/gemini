#include "CEvap.h"

CEvap* CEvap::fInstance = 0;

float const r0 = 1.16;
/**
 * Constructor
 */
CEvap::CEvap() 
{

  string fileName("tbl/evap.inp");
  string fullName;
  if (getenv("GINPUT") == NULL) fullName = fileName;
  else
    {
      string dir(getenv("GINPUT"));
     fullName = dir+fileName;
    }
  ifstream ifFile (fullName.c_str());


  if (ifFile.fail())
    {
      cout << " file " << fullName << " not found" << endl;
      abort();
    }   

  // read in number of evaporation channels 
  ifFile >> nLight;
  prob = new float [nLight];
  // create array of CLightP pointers 
  lightP = new CLightP * [nLight];

  tlArray = new CTlBarDist * [nLight]; 
  sigBarDist = new CSigBarDist * [nLight]; 

  // read in channel information and initialize pointers
  decay = new SDecay [nLight];

  // skip line
  string line;
  getline(ifFile,line);
  getline(ifFile,line);



  
  float fJ, fEx, Ek, suppress;
  int iZ, iA;
  string nameOld("");
  string name("");
  bool newTl;
  nTl = 0;
  for (int i=0;i<nLight;i++)
    {
      nameOld = name;
      ifFile >> iZ >> iA >> fJ >> fEx >> name >> suppress >> Ek;


      // determine if we need to create new transmission coeff
      newTl = 1;
      if (i > 0)
	{
	  if (name == nameOld) newTl = 0;
	}
      if (newTl)
	{
          tlArray[nTl] = new CTlBarDist(name);
          sigBarDist[nTl] = new CSigBarDist(name,(float)iZ,(float)iA);
          nTl++; 
	}       
      lightP[i] = new CLightP(iZ,iA,fJ,tlArray[nTl-1],sigBarDist[nTl-1]);
      lightP[i]->rLight = pow((float)iA,(float)(1./3.))*r0;
      lightP[i]->fEx = fEx;
      lightP[i]->suppress = suppress;
      // if an excited state add excitation energy to mass excess
      if (fEx > 0.0) lightP[i]->fExpMass += fEx;
      maxZ = iZ;
      decay[i].Ek = Ek;
      if (Ek == 0.) continue;
      // read in decay information 
        ifFile >> decay[i].Z1 >> decay[i].A1 >> decay[i].S1 >> 
	  decay[i].S2 >> decay[i].L >> decay[i].lPlusS1 >> decay[i].gamma;

	// here are sum checks to make sure decay information is possible
        if (decay[i].lPlusS1 > decay[i].S1 + (float)decay[i].L)
	  {
	  cout << "bad LPlusS1 in evap.cpp for mode " << i << endl;
	  abort();
	  }

        if (decay[i].lPlusS1 < fabs(decay[i].S1 - (float)decay[i].L))
	  {
	  cout << "bad LPlusS1 in evap.cpp for mode " << i << endl;
	  abort();
	  }

        if (fJ > decay[i].lPlusS1+decay[i].S2)
	  {
	  cout << "bad fJ in evap.cpp for mode " << i << endl;
	  abort();
	  }

        if (fJ < fabs(decay[i].lPlusS1-decay[i].S2))
	  {
	  cout << "bad fJ in evap.cpp for mode " << i << endl;
	  abort();
	  }
  
    }

  ifFile.clear();
  ifFile.close();

  
}

CEvap* CEvap::instance()
{
    if (fInstance == 0) {
        fInstance = new CEvap;
    }
    return fInstance;
}

//**************************************************
/**
 * Destructor
 */
CEvap::~CEvap()
{
  for (int i=0;i<nLight;i++) delete lightP[i];
  for (int i=0;i<nTl;i++) 
    {
     delete tlArray[i];
     delete sigBarDist[i];
    }
  delete [] prob;
  delete [] lightP;
  delete [] tlArray;
  delete [] sigBarDist;
  delete [] decay;  
}
