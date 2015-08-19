#ifndef yrast_
#define yrast_
#include <iostream>
#include <fstream>
#include "CMass.h"
using namespace std;

/**
 *!\brief fission barriers, rotational energies, etc
 *
 * Class to return barriers and rotational energies. It is constructed 
 * from bits of code such as Arnie Sierk Barfit - which gives his Finite-range
 *fission barrier, rotational energies, moments of inertia, surface areas.
 * Also I have used code from the RLDM model. In addition I enclude a 
 * interpolation of conditional fission barriers, using barriers calculated
 * by Sierk
 */

class CYrast
{
 private:

  CYrast();
  static CYrast *fInstance; //!< instance member to make this class a singleton
  static double const pi; //!< 3.14159
  //needed by getYrastRLDM
  static float const x1h[11][6]; //!< number for RLDM
  static float const x2h[11][6]; //!< numbers for RLDM
  static float const x3h[20][10]; //!< number for RLDM
  static float const x1b[11][6]; //!<numbers for RLDM
  static float const x2b[11][6]; //!<numbers for RLDM
  static float const x3b[20][10]; //!<numbers for RLDM
  //needed by Sierk functions
  static double const emncof[4][5]; //!< used in Sierk functions
  static double const elmcof[4][5]; //!<used in Sierk functions
  static double const emxcof[5][7];//!<used in Sierk functions
  static double const elzcof[7][7];//!<used in Sierk functions
  static double const egscof[5][7][5];//!<used in Sierk functions
  static double const aizroc[5][6];//!<used in Sierk functions
  static double const ai70c[5][6];//!<used in Sierk functions
  static double const ai95c[5][6];//!<used in Sierk functions
  static double const aimaxc[5][6];//!<used in Sierk functions
  static double const ai952c[5][6];//!<used in Sierk functions
  static double const aimax2c[5][6];//!<used in Sierk functions
  static double const aimax3c[4][4];//!<used in Sierk functions
  static double const aimax4c[4][4];//!<used in Sierk functions
  static double const bizroc[4][6];//!<used in Sierk functions
  static double const bi70c[4][6];//!<used in Sierk functions
  static double const bi95c[4][6];//!<used in Sierk functions
  static double const bimaxc[4][6];//!<used in Sierk functions
  static double const b[8][5][5];//!<used in Sierk functions
  void lpoly(double,int,double*);

  double A; //!< mass number
  double Z; //!< proton number
  double zz; //!< used in Sierk functions
  double amin; //!< lower limits of application of Sierk routine
  double amax; //!< upper limits of application of Sierk routine
  double pa[7]; //!<used in Sierk routines
  double pz[7];  //!<used in Sierk routines
  //needed by saddlefit
  float c[6][8][2][11][2]; //!< coeff for sadfits
  float cubic(float,float,float,float,float,float);

  static bool first; //!< only write out barrier warning once
  int Narray; //!< number of elements in array of asymmetric barriers
  static float const hbarc; //!< used for asymmetric barriers
  static float const alfinv;//!< used for asymmetric barriers
  static float const srznw;//!< used for asymmetric barriers
  static float const aknw;//!< used for asymmetric barriers
  static float const bb;//!< used for asymmetric barriers
  static float const um;//!< used for asymmetric barriers
  static float const elm;//!< used for asymmetric barriers
  static float const spdlt;//!< used for asymmetric barriers
  static float const asnw;//!< used for asymmetric barriers
  static float const kx[8];//!< used for asymmetric barriers
  static float const ky[6];//!< used for asymmetric barriers
  static float const ka[11];//!< used for asymmetric barriers
  static float const r0; //!< radius parameter
  static float const sep; //!<separation id fm between fragments
  static bool bForceSierk; //!<separation id fm between fragments
  static double addBar; //!<extrapolated Sierk barrier increase by this amount

  float sadArray[300]; //!< array stores the conditional saddle energies
  float sadArrayZA[300]; //!< array stores saddle energies after correction
  CMass * mass; //!< class for mass defects
  static float const deltaJ; //!< used to extend sierk barrier to higher J
  static float const kRotate; //!< constant for rotional energy
  int iZ; //!< proton number
  int iA; //!<mass number
  float fJ; //!< spin
 public:
  static CYrast *instance(); //!< instance member to make this class a singleton
  double Jmax;  //!< max spin where the fission barrier exists
  float getYrast(int,int,float);
  float getYrastModel(int,int,float);
  float getYrastRLDM(int,int,float);
  float getYrastSierk(float);
  float getJmaxSierk(int,int);
  float getBarrierFissionSierk(float);
  float getSymmetricSaddleEnergy(int,int,float);
  float getBarrierFissionRLDM(int,int,float);
  float getBsSierk(float);
  static void forceSierk(bool=1);
  static void printParameters();
  void prepareAsyBarrier(int, int, float);
  void printAsyBarrier();
  float getSaddlePointEnergy(int,int);
  float getSaddlePointEnergy(float);
  float getMomentOfInertiaSierk(float);
  float WignerEnergy(int iZ, int iA);

  float momInertiaMin; //!< minimum saddle-point moment of inertia 
  float momInertiaMid; //!< intermediate saddle-point moment of inertia
  float momInertiaMax; //!< maximum saddle-point moment of inertia
  // float sumPair;
  //float sumShell;
  //float viola(float,float,float,float);
};
#endif
