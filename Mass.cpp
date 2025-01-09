#include "CMass.h"
#include <cmath>

CMass* CMass::fInstance = 0;  // mod-TU


/**
 * Constructor
 */
CMass::CMass()
{
  chart = CChart::instance();
  fExpMass = new float [chart->iMassDim];
  fCalMass = new float [chart->iMassDim];
  fFRM = new float [chart->iMassDim];
  fPair = new float [chart->iMassDim];
  //fBeta2 = new float [chart->iMassDim];
  fShell = new float [chart->iMassDim];
  fShell2 = new float [chart->iMassDim];
  ReadThomasFermiFile();
  ReadFRDMFile();
  fExpMass[0] = 8.071;
  fExpMass[1] = 7.289;
  fCalMass[0] = 8.071;
  fCalMass[1] = 7.289;


  int iZmin = 4;
  int iZmax = 136;
  for (int iZ =iZmin;iZ<=iZmax;iZ++)
    {
      int iNmin = chart->getAmin(iZ)-iZ;
      int iNmax = chart->getAmax(iZ)-iZ;
      //cout << iZ << " " << iNmin << " " << iNmax << endl;
      for (int iN = iNmin;iN<= iNmax;iN++)
	{
          
          int index = chart->getIndex(iZ,iZ+iN);
          fPair[index] = getPairing2(iZ,iZ+iN);

          //redefine the shell correction to be difference between 
	  //experimental mass and (FRM masses plus given pairing corrections)
          fShell2[index] = fExpMass[index] - fPair[index] - fFRM[index];

	}
    }



  for (int iZ =iZmin;iZ<=iZmax;iZ++)
    {
      int iNmin = chart->getAmin(iZ)-iZ;
      int iNmax = chart->getAmax(iZ)-iZ;

      for (int iN = iNmin;iN<= iNmax;iN++)
	{
          int index = chart->getIndex(iZ,iZ+iN);

          float fpairN;
          if(iN%2 == 1) fpairN = 0.;
          else if (iN-1 >= iNmin && iN+1 <= iNmax)
	    {
             fpairN = -(getShellCorrection2(iZ,iZ+iN-1) + 
                        getShellCorrection2(iZ,iZ+iN+1))/2.
	    + getShellCorrection2(iZ,iZ+iN);
	    }
          else if (iN-1 < iNmin) fpairN = -getShellCorrection2(iZ,iZ+iN+1)
	    + getShellCorrection2(iZ,iZ+iN);
          else if (iN+1 > iNmax) fpairN = -getShellCorrection2(iZ,iZ+iN-1)
	    + getShellCorrection2(iZ,iZ+iN);


          //then proton pairing
          float fpairZ;
          if(iZ%2 == 1) fpairZ = 0.;
          else if (iZ+1 <= iZmax) fpairZ = 
           -(getShellCorrection2(iZ-1,iZ-1+iN) 
           + getShellCorrection2(iZ+1,iZ+iN+1))/2.
	  + getShellCorrection2(iZ,iZ+iN);
	  else if (iZ+1 > iZmax)  fpairZ = 
           -getShellCorrection2(iZ+1,iZ+iN+1)+ getShellCorrection2(iZ,iZ+iN);

	  //fpairZ = 0.;


	  fPair[index] = fpairN + fpairZ + getPairing2(iZ,iZ+iN);

          fShell[index] = fExpMass[index] - fFRM[index]  - fPair[index];


	}
    }
  


}

CMass* CMass::instance() // mod-TU
{
    if (fInstance == 0) {
        fInstance = new CMass;
    }
    return fInstance;
}

//**********************************************
/**
 * Destructor
 */
CMass::~CMass()
{
  delete [] fExpMass;
  delete [] fCalMass;
  delete [] fFRM;
  // delete [] fBeta2;
  delete [] fShell;
  delete fShell2;
}
//********************************************
/**
 * Returns the experimental mass excess
 *
 * If teh experimental excess is not known, then the Moller Nix value
 * is returned
 \param iZ is the proton number
 \param iA is the mass number
 */
float CMass::getExpMass(int iZ, int iA)
{
 
  //find location of nuclide in array where the mass is stored
  int i = chart->getIndex(iZ,iA);
  return fExpMass[i];

}
//********************************************
/**
 * Returns the calculated mass excess from Moller and Nix
 *
 \param iZ is the proton number
 \param iA is the mass number
 */
float CMass::getCalMass(int iZ, int iA)
{
 
  //find location of nuclide in array where the mass is stored
  int i = chart->getIndex(iZ,iA);
  return fCalMass[i];

}
//********************************************
/**
 * Returns the shell correction from Moller and Nix
 *
 \param iZ is the proton number
 \param iA is the mass number
 */
float CMass::getShellCorrection(int iZ, int iA)
{
 
  //find location of nuclide in array where the mass is stored
  int i = chart->getIndex(iZ,iA);
  return fShell[i];

}

//********************************************
/**
 * Returns the shell correction from Moller and Nix
 *
 \param iZ is the proton number
 \param iA is the mass number
 */
float CMass::getShellCorrection2(int iZ, int iA)
{
 
  //find location of nuclide in array where the mass is stored
  int i = chart->getIndex(iZ,iA);
  return fShell2[i];

}
//********************************************
/**
 * Returns the liquid drop mass from moller and Nix
 \param iZ is the protom number
 \param iA is the mass number
 */
float CMass::getLiquidDropMass(int iZ, int iA)
{
 
  //find location of nuclide in array where the mass is stored
  int i = chart->getIndex(iZ,iA);
  if ( i == -1) 
    { 
      return -1000;
    }
  return fFRM[i];

}
//*************************************************
  /**
   *
   * Calculates macroscopic finite range model masses of spherical
   * nucleus using formula of Krappe, Nix, and Sierk.
   *
   * Reference- (Phys Rev C20, 992 (1979))
   * modified to use the parameters of Moller + Nix Nucl. Phys. A361(1981)
   * 117. Pairing correction term for odd-odd nuclei
   * is included, as this is the most appropriate ground state for hot nuclei
   * where shell and pairing effects have washed out.
   \param iZ is the proton number
   \param iA is the mass number
   */
 float CMass::getFiniteRangeMass(int iZ, int iA)
{
  return getFiniteRangeMass((float)iZ,(float)iA);
}
//**********************************************************
  /**
   *
   * Calculates macroscopic finite range model masses of spherical
   * nucleus using formula of Krappe, Nix, and Sierk.
   *
   * Reference- (Phys Rev C20, 992 (1979))
   * modified to use the parameters of Moller + Nix Nucl. Phys. A361(1981)
   * 117. Pairing correction term for odd-odd nuclei
   * is included, as this is the most appropriate ground state for hot nuclei
   * where shell and pairing effects have washed out.
   \param fZ is the proton number
   \param fA is the mass number
   */
float CMass::getFiniteRangeMass(float fZ, float fA)
{
  float fN = fA - fZ;
  float fA13 = pow(fA,(float)(1./3.));
  float fA23 = pow(fA13,2);

  // relative neutron excess
  float fI = (fN-fZ)/fA;
  float fFiss = pow(fZ,2)/fA13;

  // neutron-proton terms
  float const fMassN = 8.071431;
  float const fMassH = 7.289034;
  float fEnz = fMassN*fN + fMassH*fZ;

  //Volume energy 
  float const fAv = 15.9937;
  float const fKv = 1.927;
  float fEvol = -fAv*(1.-fKv*pow(fI,2))*fA;

  // Surface energy
  float const fa = 0.68;
  float const fR0 = 1.16;
  float const fAs = 21.13;
  float const fKs = 2.3;
  float fX=fa/(fR0*fA13);
  float fact=1.-3.*pow(fX,2)
    + (1.+1./fX)*(2.+3.*fX+3.*pow(fX,2))*exp(-2./fX);
  float fEsurf = fAs*(1.-fKs*pow(fI,2))*fA23*fact;

  //Coulomb energy
  float const e2 = 1.4399764;
  fact = fFiss-0.76361*pow(fZ,(float)(4./3.))/fA13;
  float fECoul = 0.6*e2/fR0*fact;

  //Wigner term
  float const fW = 36.;
  float const ael = 1.433e-5;
  float fEwigner = fW*fabs(fI)-ael*pow(fZ,(float)2.39);

  //correction to Coulomb energy for diffuse surface
  //see Davies & Nix Phys. Rev. C14 (1976) 1977
  float const b = 0.99;
  float ad=0.7071*b;
  fX = ad/(fR0*fA13);
  fact= 1 - 1.875*fX+2.625*pow(fX,3) 
    -.75*exp(-2./fX)*(1.+4.5*fX+7.*pow(fX,2)+3.5*pow(fX,3));
  float fEcd = -3.*pow(fZ,2)*e2*pow(ad,2)/pow(fR0*fA13,3)*fact;

  // correction to coulomb energy due to proton form factor
  float const rp = 0.8;
  float akf=1./fR0*pow(7.06858*fZ/fA,1./3.);
  fX = pow(rp*akf,2);
  float fEcpf=-0.125*pow(rp,2)*e2/pow(fR0,3) 
    *(3.0208-0.113541667*fX+0.0012624*pow(fX,2))*pow(fZ,2)/fA;

  //A0 term
  float const c0 = 4.4;
  float fEa0=c0;

  //Charge asymmetry term
  float const ca = 0.212;
  float fEca=ca*(pow(fZ,2)-pow(fN,2))/fA;

  // pairing term for odd-odd nuclei
  float deltau = 12./sqrt(fA);
  float deltal = 20./fA;
  float fEpair =deltau - 0.5 * deltal;

  // add all terms
  return fEnz+fEvol+fEsurf+fECoul+fEwigner+fEcd+fEcpf+fEa0+fEca+fEpair;
}
//******************************************
/**
 * Returns the pairing correction to the mass formula.
 * From from Moller Nix is used.
 \param iZ is the proton number
 \param iA is the mass number
 */
float CMass::getPairing(int iZ, int iA)
{

  return fPair[chart->getIndex(iZ,iA)];
}

//****************************************
/**
 * Returns the pairing correction to the mass formula.
 * From from Moller Nix is used.
 \param iZ is the proton number
 \param iA is the mass number
 */
float CMass::getPairing2(int iZ, int iA)
{


  if(iZ==0 || iZ==iA)
    return 0.;
  float fZ = iZ;
  float fA = iA;
  float fN = fA - fZ;
  int iN = iA - iZ;
  int  ioez = iZ%2;
  int ioen = iN%2;
  float fPairing;
  if (iN == iZ && ioez ==1)
    {
      fPairing = 4.8/pow(fN,(float)(1./3.)) 
	+ 4.8/pow(fZ,(float)(1./3.)) - 6.6/pow(fA,(float)(2./3.)) 
       + 30./fA;
    }
  else if (ioez == 1 && ioen == 1)
    {
      fPairing = 4.8/pow(fN,(float)(1./3.)) 
	+ 4.8/pow(fZ,(float)(1./3.)) - 6.6/pow(fA,(float)(2./3.));     
    }
  else if (ioez == 1 && ioen == 0)
    {
      fPairing = 4.8/pow(fZ,(float)(1./3.));
    }
  else if (ioez == 0 && ioen == 1)
    {
      fPairing = 4.8/pow(fN,(float)(1./3.));
    }
  else fPairing = 0.;

  //want to redefine odd-odd to have zero paring energy
  fPairing -= 4.8/pow(fN,(float)(1./3.)) 
	+ 4.8/pow(fZ,(float)(1./3.)) - 6.6/pow(fA,(float)(2./3.));     

  if (iN == iZ) fPairing -= 30./fA;

  return fPairing;

  
}
//******************************************
/**
 * Reads in the mass table from Moller and Nix
 */
void CMass::ReadFRDMFile()
{


  string fileName("tbl/mass.tbl");
  string fullName=string(GINPUT)+fileName;
  ifstream ifFile (fullName.c_str());


  if (ifFile.fail()) 
    {
      cout << "could not open" << fullName << endl;
      abort();
    }


  while (!ifFile.eof())
    {
      char integ[6]={"     "};
      int izz,iaa;
      for (int i=0;i<5;i++) integ[i] = ifFile.get();
      izz = atoi(integ);
      for (int i=0;i<5;i++) integ[i] = ifFile.get();
      // inn = atoi(integ); // Unused

      for (int i=0;i<5;i++) integ[i] = ifFile.get();
      iaa = atoi(integ);

      char floaty[11]={"          "};
      for (int i=0;i<10;i++) floaty[i] = ifFile.get();
      for (int i=0;i<10;i++) floaty[i] = ifFile.get();
      for (int i=0;i<10;i++) floaty[i] = ifFile.get();
      // fb = atof(floaty); // Unused
      for (int i=0;i<10;i++) floaty[i] = ifFile.get();
      for (int i=0;i<10;i++) floaty[i] = ifFile.get();
      for (int i=0;i<10;i++) floaty[i] = ifFile.get();
      for (int i=0;i<10;i++) floaty[i] = ifFile.get();
      for (int i=0;i<10;i++) floaty[i] = ifFile.get();
      for (int i=0;i<10;i++) floaty[i] = ifFile.get();

      float f1,f2,f3;
      for (int i=0;i<10;i++) floaty[i] = ifFile.get();
      f1 = atof(floaty);
      for (int i=0;i<10;i++) floaty[i] = ifFile.get();
      f2 = atof(floaty); //calculated mass
      bool there = 0;
      for (int i=0;i<10;i++) 
	{
         floaty[i] = ifFile.get();
	 if (floaty[i] != ' ') there = 1;
	}
      f3 = atof(floaty);//experimental mass
      for (int i=0;i<10;i++) floaty[i] = ifFile.get();
      // f4 = atof(floaty); // Unused
      for (int i=0;i<10;i++) floaty[i] = ifFile.get();
      // f5 = atof(floaty); // Unused
      for (int i=0;i<10;i++) floaty[i] = ifFile.get();
      // f6 = atof(floaty); // Unused


      //read to end of line
      for (;;) 
	{
	  if (ifFile.get() == 10) break;
          if (ifFile.eof()) break;
	}
      int index = chart->getIndex(izz,iaa);
      if (index >= 0)
	{
          float fPair = getPairing2(izz,iaa);

	  if (there) fExpMass[index] = f3;
	  else fExpMass[index] = f2;
          fCalMass[index] = f2;
          fFRM[index] = f2 - f1 - fPair;
          fShell[index] = f1;
          //fBeta2[index] = fb; 


	  /*	  
	  if (izz == 10) cout << iaa << " " << iaa-izz << " " << 
	    fFRM[index] << " " << fShell[index] << " " << fPair <<
            " " << fExpMass[index] << " " << fCalMass[index] << " " 
			      << there << endl;
	  */

          
	}

    }

   ifFile.clear();
   ifFile.close();
}

//******************************************
/**
 * Reads in the mass table from the Thomas Fermi Model
 * of Myers and Swietcki
 */
void CMass::ReadThomasFermiFile()
{


  string fileName("tbl/mass_tf.tbl");
  string fullName=string(GINPUT)+fileName;
  ifstream ifFile (fullName.c_str());

  if (ifFile.fail()) 
    {
      cout << "could not open mass_tf.tbl" << endl;
      abort();
    }

  //skip  introductory remarks in file
  string line;
  for (int i=0;i<66;i++) 
    {
      getline(ifFile,line);
    }
  while (!ifFile.eof())
    {
      char integ3[4]={"   "};
      
      int izz,iaa;
      for (int i=0;i<3;i++) integ3[i] = ifFile.get();
      char integ4[5]={"    "};
      izz = atoi(integ3);

      for (int i=0;i<4;i++) integ4[i] = ifFile.get();
      //int inn = atoi(integ4); //Unused
      for (int i=0;i<4;i++) integ4[i] = ifFile.get();
      iaa = atoi(integ4);



      char float5[11]={"          "};
      for (int i=0;i<5;i++) float5[i] = ifFile.get();

      for (int i=0;i<5;i++) float5[i] = ifFile.get();
      // char float8[9]={"        "}; // Unused
      for (int i=0;i<8;i++) float5[i] = ifFile.get();

      char float6[7]={"      "};
      for (int i=0;i<6;i++) float6[i] = ifFile.get();
      float fshell = atof(float6);

      for (int i=0;i<5;i++) float5[i] = ifFile.get();
      float fPairing = atof(float5);
      for (int i=0;i<6;i++) float6[i] = ifFile.get();
      for (int i=0;i<8;i++) float5[i] = ifFile.get();
      char float7[8]={"       "};
      for (int i=0;i<7;i++) float7[i] = ifFile.get();
      float f1 = atof(float7);
      for (int i=0;i<7;i++) float7[i] = ifFile.get();
      float f2 = atof(float7);

      //read to end of line
      for (;;) 
	{
	  if (ifFile.get() == 10) break;
          if (ifFile.eof()) break;
	}

      int index = chart->getIndex(izz,iaa);
      if (index >= 0)
	{
          //note do not use fPairing from table,
	  //we have refined odd-odd to have zero pairing
          fPairing = getPairing2(izz,iaa);
          fExpMass[index] = f2;
          if (izz < 5) fCalMass[index] = f2;
          else fCalMass[index] = f1;
          fFRM[index] = fCalMass[index] - fshell - fPairing;
          fShell[index] = fshell*2.;

	}
      
    }

   ifFile.clear();
   ifFile.close();
}

