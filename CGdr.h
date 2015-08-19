#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;

/**
 *!\brief stores parameters for a component of the GDR lineshape
 *
 */
struct Lshape
{
  float strength; //!< Lorentzian strength (fraction)
  float energy;   //!< Centroid [MeV]
  float gamma;    //!< Width [MeV]
};


/**
 *!\brief user-defined GDR line shape
 *
 * Reads in and calculates the GDR line shape in the E1 gamma-decay
 * branch of a hot compound nucleus. It reads file tbl/GDR.inp 
 * for the parameters of up to 5 Lorentzians.
 * The parameters are the strength, the centriod [MeV] and the 
 * width [MeV]. The sum of the stregths should add to unity 
 * 
 */

class CGdr
{
 public: 
  static CGdr* instance(); //!< instance member to make this a singleton
  float getLineShape(float e);

 protected:
  static  CGdr* fInstance; //!< instance member to make this a singleton
    CGdr();                //!< constructor
    int N;                 //!< number of Lorentzians in lineshape
    Lshape lineShape[5];   //!< parameters for each Lorentzian
};
