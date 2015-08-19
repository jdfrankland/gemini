#include <fstream>
#include <iostream>
#include <cstdlib>


using namespace std;
/**
 *!\brief storage of chart of isotope
 *
 * contain max and min A values for a given Z which are allowed in GEMINI
 */
struct SIsotope
{
  short unsigned iAmax; //!< maximun A value
  short unsigned iAmin; //!< minimum A value
};



/**
 *!\brief defines accessable region of chart of nuclides

 * contains the neutron rich and proton rich 
 * limits to the chart of nuclides for which 
 * decay fragments can be chosen
 */

class CChart
{
 private:
  CChart();
  static CChart *fInstance; //!< instance member to make this class a singleton
  SIsotope * isotope; //!< lists max and min A for each Z
  static int const iZmax; //!< max Z allowed in GEMINI
  int * iZindex; //!< list the array number of the lightest isotope of each Z  
  

 public:
 static CChart* instance(); //!< instance member to make this a singleton
 ~CChart();
 int getAmax(int iZ);
 int getAmin(int iZ);
 int getIndex(int iZ,int iA);
 int iMassDim; //!< dimension of mass array
};
