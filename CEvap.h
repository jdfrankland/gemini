#ifndef evap_
#define evap_

#include "CLightP.h" 
#include <fstream>
#include <iostream>

using namespace std;


/**
 *!\brief storage of light particle properties
 *
 * the Structure stores the information on the secondary decay of
 * particle-unstable evaporation fragments. (used in CEvap) 
 */
struct SDecay
{
  float Ek; //!< kinetic energy released in secoardy decay (MeV)
  float S1; //!< spin of one of the secondary particles
  float S2; //!< spin of the other secondary particle
  float lPlusS1; //!< spin + orbital AM of the first secondary particle
  short unsigned Z1; //!<proton number of the first secondary particle
  short unsigned A1; //!< mass number of the first secondary particle
  short unsigned L; //!< orbital angular momentum of the decay
  float gamma; //!< width of the state im MeV
};


/**
 *!\brief stores info for light particle evaporation
 *
 * class to store indormation on all the light particle evaporation channels
 */
class CEvap
{
 protected:
  CEvap();
  static CEvap *fInstance; //!< instance member to make this class a singleton
  CTlBarDist ** tlArray; //!< arrays of objects generating transmission coef.
  CSigBarDist ** sigBarDist; //!< arrays of objects generating inverse cross section
  int nTl; //!< number of transmission coef generating objects

 public:
  static CEvap *instance();
  int nLight; //!< number of light-particle evaporation channels
  CLightP ** lightP; //!< array of evaporated particles
  int maxZ; //!< max Z of all fragments evaporated, IMF emission for larger Z's
  SDecay * decay; //!<information on secondary decay of evaporated particles

  ~CEvap();
  float * prob; //!< probability for decay
};

#endif
