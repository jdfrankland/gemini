// -*- mode: c++ -*- 
//
#ifndef nucleus_
#define nucleus_

#include "CMass.h"
#include "CYrast.h"
#include "CLevelDensity.h"
#include "CAngle.h"
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "CAngle.h"
#include "CEvap.h"
#include "CAngleDist.h"
#include "CScission.h"
#include "CWeight.h"
#include "SStoreEvap.h"
#include <vector>
#include "CGdr.h"
using namespace std;

/**
 *!\brief storage
 *
 * this structure store information of the evaporated particle orbital AM
 */
struct SStoreSub
{
  float gamma; //!< weight factor for the orbital amgular momentum
  short unsigned L; //!< orbital angular momentum
};

typedef vector<SStoreSub> SStoreSubVector;
typedef vector<SStoreSub>::const_iterator SStoreSubIter;

/**
 *!\brief storage
 *
 * this structure info on complex fragment channels
 */

struct SStore
{
  double gamma; //!< decay width of channel
  short unsigned iZ; //!< proton number of complex fragment
  short unsigned iA; //!< mass number of complex fragment
};



typedef vector<SStore> SStoreVector;
typedef vector<SStore>::const_iterator SStoreIter;

/**
 *!\brief storage
 *
 * structure stores information of evaporation decay sub channels
 */
struct SStoreChan
{
  float S2; //!< spin of daughter
  float Ek; //!< kinetic energy of evaporated particle (MeV)
  float Ex; //!< excitation energy of daughter (MeV)
  float gamma; //!< partial decay width for this subchannel (MeV)
  float temp; //!< temperature of daughter (MeV)
  float UPA; //!< thermal excitation energy per nucleus of daughter (MeV)
  short unsigned L; //!< orbital angular momentum of evaporated particle
};


/**
 *!\brief Hauser-Feshbach, Bohr-Wheeler, Morretto, formulisms
 *
 * Class CNucleus impliments the GEMINI statistical mode code.
 * It follows the decay of a compound nucleus by a sequential series of 
 * binary decays. The decay widths are calculated with Hauser-Feshbash
 * formulism for light particles and Morreto's transition-state formulism
 * for other binary decays
 */


class CNucleus : public CNuclide, public CWeight
{
 protected:
  float Ecoul; //!< Coulomb barrier (HauserFeshbach)

  bool notStatistical; //!< this does not decay statistically, evap. frag. only

  short unsigned notStatisticalMode;//!< specifies type of nonStatisical decay
  float fPairing; //!< pairing energy
  float fShell; //!<shell correction
  float fU0; //!< thermal excitation energy
  float Erot; //!<yrast energy
  float Jmax; //!< max spin with a fission barrier
  float fMInertia; //!< spherical moment of inertia
  float logLevelDensity; //!< store the log of the level density of the nucleus
  float temp; //!< nuclear temperature
  int fissionZ; //!< proton number of fission fragment
  int fissionA; //!< mass number of fission fragment
  int fissioningZ; //!< proton number of fission parent
  int fissioningA; //!< mass number of fission parent
  int iZ1_IMF_Max; //!< maximum Z for IMF emission

  float fissionU; //!< thermal excitation energy of both fission fragments
  float EdefScission; //!< deformation energy of the scission configuration

  bool saddleToSciss; //!< indicated decay during saddle-to-scission transition
  float timeSinceSaddle; //!< stores the time since the saddle was crossed
  float timeSinceStart; //!< stores the time since the decay began

  void saddleToScission();
  void massAsymmetry(bool);
  bool needSymmetricFission; //!< indicated the Bohr-Wheeler width is needed


  static bool const noSymmetry;//!< true - old gemini with Morreto for all 
  float timeScission; //!< time required to go from saddle to scission
  static float const viscosity_scission; //!< viscosity during saddleTosciss
  static float const viscosity_saddle; //!< viscosity  during saddleTosciss
  static float timeTransient; //!< transient fission delay 
  static float fissionScaleFactor; //!< fission width scaled by this factor
  static float barAdd; //!< adds to Sierk fission barrier
  static unsigned iPoint; //!< pointer to array of stable fragments
  static int iHF; //!< set evaporation mode 
  int HF; //!< evaporation mode chosen for a given decay
  static bool noIMF; //!< no imf emission is considered
  static bool BohrWheeler; //!< no imf emission is considered
  float selectJ(float,float,float,float);

  static short unsigned Zshell; //!< enforce shell effects in evaporation
  static CYrast *yrast; //!< gives fission barriers and rotational energies
  CScission scission; //!< gives scission energeis, etc
  static CLevelDensity *levelDensity; //!< gives level densities
  static CGdr * GDR ; //!< uder defined GDR line shape
  bool  bStable; //!< indicated this nucleus is particle-stable
  static float const r0; //!< radius const (fm)
  static float const sep; //!< separation between fragments
  static float threshold; //!< used to turn off unlikey evaporations

  CAngle spin; //!< orientation of the spin axis
  float velocity[3]; //!< velocity vector of nucleus in cm/ns 
  float momentum[3]; //!< momentum vector in MeV/c


  static CEvap *evap; //!< stores info on evaporated particles


  CLightP * lightP; //!< points to the light-particle decay mode
  float S2Loop(float Ekvalue);
  float S2Width(float Ekvalue);
  float EkWidth(float ek);
  void getSpin( bool saddle);

  float EkLoop();
  float getSumTl(float,float);
  float getWidthZA(float,short);
  void angleEvap();
  void angleIsotropic();
  void angleGamma();
  float S2Start; //!< Hauser-Feshback spin of daughter
  float UMin; //!< min thermal excitation energy in Hauser-Feshbach
  float EcostMin; //!<the min of the energetic cost of emitting light particles
  static short unsigned const lMaxQuantum; //!< number of l-waves to store angular dist
  static float de;//!< kinetic-energy interval for integrating in Hauser-Feshb
  int lMin; //!< minimum orbital AM for Hauser-Feshbach
  int lMax;//!< maximum orbital AM for Hauser-Feshbach
  float lPlusSMax; //!< max value of l+S of evaporated particle
  float lPlusSMin; //!< min value of l+S of evaporated particle
  float rResidue; //!< radius of daughter
  float rLightP; //!< radius of evaporated particle
  //float fMInertiaOrbit; //!< moment of Inertia for orbital motion
  float S2; //!< spin of daughter
  float EYrast2; //!< rotational energy of daughter
  SStoreEvap * storeEvap; //!< information of evap sub channels
  SStoreSub * storeSub; //!< store info on l distribution

  int iStore; //!< actual number of evap sub channels
  short unsigned EvapZ2; //!< proton number of daughter after evap.
  short unsigned EvapA2; //!< mass number of daughter after evap.
  short unsigned EvapZ1; //!< proton number of evaporated particle
  short unsigned EvapA1; //!< mass number of evaporated particle
  short unsigned EvapL; //!< orbital AM of evaporated particle
  short unsigned EvapMode; //!< ID number of evap channel 
  float EvapEx1; //!< excitation ennergy of evap. particle
  float EvapEx2; //!< excitation energy of daughter after evap.
  float EvapS2; //!< spin of daughter after evap
  float EvapS1; //!< spin of evaporated particle
  float EvapEk; //!< kinetic energy of evaporated particle (MeV)
  float EvapLPlusS; //!< toatl spin plus orbital AM of evaporated particle

  static CAngleDist angleDist; //!< selects angular distributions of decays

  float GammaEx;//!< excitation energy after gamma emission
  float GammaJ; //!< spin after gamma emission
  int GammaL; //!< gamma type E1=1, E2 = 2

  static float const gammaInhibition[3]; 
//!<scaling of gamma width from Weisskopf value
  static float const wue[3]; //!<coeff for Weisskopf units (gamma decay)
  void binaryDecay();
  void exciteScission(float,float,bool sym=1);
  float asyFissionWidth();
  float asyFissionWidthZA();
  float asyFissionWidthBW();
  void force8Be();
  void force5Li();
  void force5He();
  void force9B();

  float evaporationWidthSS();
  float gammaWidth();
  float gammaWidthE1GDR();
  float gammaWidthMultipole(int);
  float hauserFeshbach(int);
  float weiskopf( bool saddle);
  void  asyFissionDivide();
  void recursiveDecay();
  CNucleus*daughterLight;  //!< pointer to the lighter of the decay products
  CNucleus*daughterHeavy;  //!< pointer to the heavier of the decay products
  CNucleus*parent; //!< pointer to the parent nucleus

  static int const Nproducts; 
  //!< total number of possible  decay products from all decays


  static vector<CNucleus *> allProducts;
  //!< array of pointer to all decay products (stable or intermediate) 

  static vector<CNucleus *> stableProducts;
  //!< array of pointers to all stable decay products for all CN decays


  bool bResidue; //!< true if decay produced an evaporation residue
  bool bSymmetricFission; //!< true if decay resulted in symmetric fission
  bool bAsymmetricFission; //!< true if decay resulted in asymmetric fission
  int multPostLight; //!< number of post-fission neutrons for lighter ff
  int multPostHeavy; //!< number of post-fission neutrons for heavier ff
  int multPreSaddle; //!< number of pre-scission neutrons emitted
  int multSaddleToScission; 
  //!< number of neutrons emitted between saddle and scission


  static float const kRotate; //!< constant to calculated rotational energy
  void split(CAngle);

  float sigma2; //!< variance of fission mass distribution
  float symSaddlePoint;//!< symmetric saddle point energy
  static float sumGammaEnergy; //!< store the energy emitted in gamma rays

  static  vector <float> GammaRayEnergy; //!< store each gamma ray energy
  static int  nGammaRays;  //!< number of emitted gamma rays

  static bool  GDRParam; //!< if true, the standard formula for GDR decay width is used, if false the parametrized version




 public:



  bool abortEvent; //!< abort the event
  float evaporationWidth();
  float BohrWheelerWidth();
  float LestoneFissionWidth();
  float LestoneCorrection(float Usaddle, float momInertiaEff,short iAfAn);
  static float const pi; //!< 3.14159
  static double const EkFraction; // !< calculates the Ek spectra down to this
                                 // fraction of the maximum
  //functions
  CNucleus();
  CNucleus(int iZ,int iA);
  CNucleus(int iZ,int iA, float fEx, float fJ);
  ~CNucleus();
  float getSumGammaEnergy();
  int   getnGammaRays(); //get number of emitted gamma rays
  float getGammaRayEnergy(int number); // get gamma ray energy
  float getTime();

  void setNewIsotope(int iZ0, int iA0, float fEx0, float fJ0); 

  void setCompoundNucleus(float fEx0,float fJ0);
  void setCompoundNucleus(float fEx0,double dJ0);
  void setCompoundNucleus(double dEx0,float fJ0);
  void setCompoundNucleus(double dEx0,double dJ0);

  void setSpinAxis(CAngle angle);
  void setSpinAxisDegrees(CAngle angle);
  void setVelocityPolar(float =0.,float=0.,float=0.);
  void setVelocityCartesian(float vx=0.,float vy=0.,float vz=0.);

  void reset();
  static void resetGlobal();

  void print();
  void printStableProducts();
  void printAllProducts();
  void vCMofAllProducts();
  void energyConservation();


  CNucleus* getProducts(int=-1);
  CNucleus* getParent();
  CNucleus* getLightDaughter();
  CNucleus* getHeavyDaughter();
  CNucleus* getCompoundNucleus();


  int getNumberOfProducts();
  int getZmaxEvap();


  void excite(float,float);
  void excite(float,double);
  void excite(double,float);
  void excite(double,double);
  void excite(float);

  float getTheta();
  float getThetaDegrees();
  CAngle getAngle();
  CAngle getAngleDegrees();
  float getKE();
  float getVelocity();
  float getMomentum();
  float* getVelocityVector();
  float* getMomentumVector();

  static void setTimeTransient(float time);
  static void setFissionScaleFactor(float factor);
  static void setBarWidth(float width);
  static void setDeltaE(float de0);
  static void setThreshold(float threshold0);
  static void setAddToFisBarrier(float barAdd0);
  static void setNoIMF();
  static void setYesIMF();
  static void setLestone();
  static void setBohrWheeler();
  static void setSolution(int isol);
  static void setEvapMode(int iHF0=2);
  static void setUserGDR(bool mode = true);

  static float getTimeTransient();
  static float getFissionScaleFactor();
  static float getBarWidth();
  static float getDeltaE();
  static float getThreshold();
  static float getAddToFisBarrier();

  void decay();
  bool isAsymmetricFission();
  bool isSymmetricFission();
  bool isNotStatistical();
  bool isSaddleToScission();
  bool isResidue();
  int getMultPost();
  int getMultPre();
  int getMultPostLight();
  int getMultPostHeavy();
  int getMultPreSaddle();
  int getMultSaddleToScission();
  float getFissionTimeSymmetric(float & timeScission);
  float getFissionTimeAsymmetric();
  float getDecayWidth();
  float getLogLevelDensity();
  int origin; //!< specifies the origin of the fragment, prefission, post , etc
  int origin2; //!< specifies the origin of the fragment, prefission, post , etc
  void printParameters();


  //*******ROOT********
    //ClassDef(CNucleus,1)  //Gemini Nucleus
};
#endif 
