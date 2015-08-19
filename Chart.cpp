//
//___________________________________________________________
//
//  CChart
//





#include "CChart.h"

int const CChart::iZmax = 136;

CChart* CChart::fInstance = 0;

//****************************************************
  /**
   * Constructor reads in files with neutron and 
   * proton rick limits
   */
CChart::CChart()
{
  string fileName("tbl/chart.tbl");
  string fullName;
  if (getenv("GINPUT") == NULL) fullName = fileName;
  else
    {
      string dir(getenv("GINPUT"));
      fullName = dir+fileName;
    }
  ifstream ifFile (fullName.c_str());

  if (ifFile.is_open() != 1) 
    {
      cout << "file " << fullName   << "is not found" << endl;
      abort();
    } 

  isotope = new SIsotope[iZmax+1];
  iZindex = new int [iZmax+1];
  int iZ,iAmin,iAmax;
  for (;;)
    {
      ifFile >> iZ >> iAmin >> iAmax;
      if (ifFile.bad()) break;
      if (ifFile.eof()) break;
      if (iZ > iZmax) break;
      isotope[iZ].iAmin = iAmin;
      isotope[iZ].iAmax = iAmax;
    }
  ifFile.close();
  ifFile.clear();

  //construct index file
  iMassDim = 0;
  for ( int i=0;i<iZmax+1;i++)
    {
      iZindex[i] = iMassDim;
      iMassDim += isotope[i].iAmax - isotope[i].iAmin + 1;
    }

}

CChart* CChart::instance()
{
    if (fInstance == 0) {
        fInstance = new CChart;
    }
    return fInstance;
}

//***********************************************************
  /**
   *descructor
   */
CChart::~CChart()
{
  //descructor
  delete [] isotope;
  delete [] iZindex;
}
//**************************************************
  /**
   *returns the maxium iA for a given element that will be considered 
   *in the decay 
   \param iZ is the proton number
   */
int CChart::getAmin(int iZ)
{
  //returns the minimum A value for a given Z
  if (iZ > iZmax) 
    {
      cout << "CHart above its limits" << endl;
      cout << "iZ=" << iZ << endl;
      abort();
    }
  return isotope[iZ].iAmin;
}
//**************************************************
  /**
   * Returns the minimum iA for a given element that will be considered 
   *in the decay 
  \param iZ is the proton number
  */

int CChart::getAmax(int iZ)
{
  //returns the maximum A value for a given Z
  if (iZ > iZmax) 
    {
      cout << "CChart above its limits" << endl;
      cout << "iZ=" << iZ << endl;
      abort();
    }
  return isotope[iZ].iAmax;
}
//*******************************************************
/**
 * Returns the index number of a particular nuclide
\param iZ is the proton number
\param iA is the mass number
 */
int CChart::getIndex(int iZ, int iA)
{
  //returns the index in the mass table for a specified Z,A
  if (iZ < 0) 
    {
      cout << " Z < 0 in chart" <<endl;
      return -1;
    }
  else if (iZ > iZmax)
    {
      //cout << " Z > iZmax in CChart " << endl;
      return -1;
    }
  if (iA < isotope[iZ].iAmin || iA > isotope[iZ].iAmax)
    {
      //cout << " outside chart of nuclides defined in CChart" << endl;
      return -1;
    }
  return iZindex[iZ] + iA - isotope[iZ].iAmin; 
}
