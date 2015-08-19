#ifndef tlArray_
#define tlArray_

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;


/**
 *!\brief IWBC fit parameters
 *
 * This structure contains the parameters with were obtained from a fit
 * to the IWBC transmission coefficients 
 */
struct SAngTl
{
  float  coef[7]; //!< fit parameters
};

/**
 *!\brief IWBC fit parameters
 *
 * Stores transmission coef. fit parameters for each daughter proton number
 */
struct SZcoef
{
  SAngTl Tl[121]; //!< array of fit paramters for each daughter Z
};


/**
 *!\brief IWBC transmission coefficients
 *
 * Provides transmission coefficients and inverse cross sections
 * for light-partice evaporation
 */

class CTlArray
{
protected:
  SZcoef zcoef; //!< structure to store coefficients parametrizating Tl's
  string sName; //!< name of file containg Tl parameters
  float shift; //!< energy shift used in the parametrization
  SAngTl *trans; //!<pointer to coefficents for a select iZ

 public:
  CTlArray(string);
  float getTermInExp(int iL,float fEk);
  void  prepare(int iZ);
  float getTl(int iL ,float fEk);
  float getInverseXsec(float fEk);
  int iZMin;  //!< the minimum Z-value for coefficinets are stored
};

#endif
