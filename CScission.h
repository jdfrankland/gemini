#include <cmath>
#include "CMass.h"
/**
 *!\brief scission energies
 *
 * calculates information on scission point for saddle-to-scission
 * evaporation and for fission mass distributions
 */


class CScission
{
 protected:
  static float const e2;  //!< \f$e^{2} = 1.44 MeV fm^{-1}\f$
  static float const kRotate; //!< constant for rotational energy 
  static float const slopeViola; //!< for Viola systematics of fission KE
  static float const constViola; //!< for Viola systematics of fission KE

  float mu0; //!< reduced mass
  float Vc0; //!< Coulomb self energy of a symmetric frag
  float R1; //!< radius of  frag1 in fm
  float momInertia1; //!< moment of inertia of  frgament1
  float R2; //!< radius of  frag1 in fm
  float momInertia2; //!< moment of inertia of  frgament1
  float k1; //!< part of the rotational energy
  bool sym; //!< logical symmetic of non symmetric fission
  float Z1;//!< atomic number of lighter fragment after fission
  float Z2;//!< atomic number of heavier fragment
  float A1;//!< mass number of lighter fragment after fission
  float A2;//!< mass number of heavier fragment 

  CMass * mass; //!< gives mass excess 

 public:
  static float const r0; //!< constant R=r0*A^(1/3)
  float sep; //!< separation between the surafces of the 2 spheres
  float sep0; //!< separation determined in init()
  float sep1; //!< separation determined in sigmaFissionSystematics()
  int iA; //!< mass number of system
  int iZ; //!< proton number of system

  float A; //!< float value of iA
  float Z; //!< float value of iZ
  float fJ; //!< spin of system
  float Esymmetric; //!< energy for symmetric mass split
  float ekTot; //!< total fission kinetic energy from getFissionKineticEnergy()
  float Erotate1; //!<rotational energy of fragment1
  float Erotate2; //!<rotational energy of fragment2
  float Epair1; //!< pairing energy of frag1
  float Epair2; //!< pairing energy of frag2
  float Eshell1; //!< shell energy of frag1
  float Eshell2; //!< shell energy of frag2
  float EkCoul; //!< Colomb part of total fisison kinetic energy
  float EkRot; //!< rotational part of total fission kinetic energy

  CScission(int iZ0,int iA0,float fJ0, int iChan);
  CScission();
  void init(int iZ0,int iA0,float fJ0,int ichan,float Z1=0.,float A1=0.);
  float getSep(float EkViola);
  float getScissionEnergy();
  float getScissionEnergy(int iZ1,int iA1);
  float getFissionKineticEnergy(int iZ1, int iA1);
  float sigmaFissionSystematicsScission(int iZ, int iA, float fJ, float fUscission);
  float sigmaFissionSystematicsSaddle(int iZ, int iA, float fJ, float fUscission);
};
