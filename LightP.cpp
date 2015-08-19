#include "CLightP.h"


  //! constructor which creates a CTlBarDist object
  /**
    \param iZ0 proton number of particle
    \param iA0 mass number of particle
    \param fJ0 spin of particle
    \param name0 name of *.tl file containing transmission coefs.
   */
CLightP::CLightP(int iZ0, int iA0, float fJ0, string name0) 
   : CNuclide(iZ0,iA0,name0)
{
   fJ = fJ0;

   tlArray = new CTlBarDist(name0);
   sigBarDist = new CSigBarDist(name0,(float)iZ0,(float)iA0);   

}
//******************************************************
  //! constructor which uses an existing CTlBarDist object
/**
 * Alternative constructor which uses a previously defined
 * transmission coeff object CTlArray
    \param iZ0 proton number of particle
    \param iA0 mass number of particle
    \param fJ0 spin of particle
    \param tlArray0 is the pointer to the object for transmission coefficients
    \param sigBarDist0 is the pointer to the object for inverse xsections
 */ 
CLightP::CLightP(int iZ0, int iA0, float fJ0, CTlBarDist* tlArray0,
   CSigBarDist * sigBarDist0) 
   : CNuclide(iZ0,iA0)
{
   fJ = fJ0;
   tlArray = tlArray0;
   sigBarDist = sigBarDist0;
   mustDelete = 0;

} 

//*****************************************************************
/**
 * Destructor
 */
CLightP::~CLightP()
{
  if (mustDelete)
    {
      delete tlArray;
      delete sigBarDist;
    }
}
