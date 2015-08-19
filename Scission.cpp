#include "CScission.h"

float const CScission::slopeViola = .1189;
float const CScission::constViola = 7.3;
float const CScission::e2=1.44;
float const CScission::r0=1.2;
float const CScission::kRotate = 41.563;

/**
 * simple constructor
 */
CScission::CScission()
{
  mass = CMass::instance();//mass singleton
}

/**
 * constructor
/param iZ0 is proton number of fissioning nucleus
/param iA0 is mass number of fission nucleus
/param fJ0 is spin of fissioning nucleus
/param iChan =1 for imf, =2 symmetric fission
 */
CScission::CScission(int iZ0, int iA0, float fJ0, int iChan)
{
  init(iZ0,iA0,fJ0,iChan);
}


//*******************************************************
  /**
   * initialized the calls for a given nucleus
    /param iZ0 is the proton number of the nucleus
    /param iA0 is the mass number of the nucleus
    /param fJ0 is the spin of the nucleus
    /param iChan =1 for imf, =2 symmetric fission
   */
void CScission::init(int iZ0, int iA0, float fJ0,int iChan, 
             float Z10/*=0.*/,float A10/*0.*/)
{
  iA = iA0;
  iZ = iZ0;
  fJ = fJ0;
  A = (float)iA;
  Z = (float)iZ;

  if (A10 == 0)
    {
     sym = 1;
     A1 = A/2.;
     Z1 = Z/2.;
    }
  else 
    {
     sym = 0;
     Z1 = Z10;
     A1 = A10;
    }

  A2 = (float)iA - A1;
  Z2 = (float)iZ - Z1;

  R1 = pow(A1,(float)(1./3.))*r0;
  momInertia1 = 0.4*A1*pow(R1,2);
  R2 = pow(A2,(float)(1./3.))*r0;
  momInertia2 = 0.4*A2*pow(R2,2);


  Vc0 = Z1*Z2*e2;
  mu0 = A1*A2/(A1+A2);
  k1 = kRotate/2.*mu0*pow(fJ,2);
  

  float Z2A13;
  if (sym == 1)Z2A13 = pow((float)iZ,2)/pow((float)iA,(float)(1./3.));
  else Z2A13= 6.349*Z1*Z2/(pow(A1,float(1./3.))+pow(A2,float(1./3.)));
  float ekViola = 0.;
  if (iChan == 1)
    {
     ekViola = slopeViola*Z2A13 + constViola;
    ekViola += 0.004*pow(fJ0,2);   // angular momentum Viola energy
    }
  else 
    {
     //alternatibve from Rusanov et al
     if (Z2A13 < 900.) ekViola = 0.131*Z2A13;
     else ekViola = 0.104*Z2A13 + 24.3;
    }
 

  //in a two sphere approx, go find separation energy which gives desired
  //fission kinetic energy
  sep = getSep(ekViola); 
  sep0 = sep;
}
//******************************************************

  /**
   * approximates the scission configuration by two separated spheres.
   * returns the separation between the surface of the sphere which
   * is consistent with input fission total kinetic energy
   /param ekViola is the total fission kinetic energy
   */

float CScission::getSep(float ekViola)
{
  float R = R1+R2 + 4.;
  for(;;)
    {
      float momInertiaTot = momInertia1+momInertia2+mu0*pow(R,2);
      float ek = Vc0/R + k1*pow(R,2)/pow(momInertiaTot,2);
      if (fabs(ek-ekViola) < 0.1) break;
      float dek = -Vc0/pow(R,2) - 4.*mu0*k1*pow(R,3)/pow(momInertiaTot,3)
	+ 2*R/pow(momInertiaTot,2);
      float dR = -(ek-ekViola)/dek ;
      if (R+dR < 0.) R=0.9*R;
      else R += dR;
    }
  return R - R1 - R2;
}
//***************************************************
  /**
   * returns the symmetric fission scission energy
   * from a two sphere approximation using the separation sep
   * previously determined from init() to sigmaFissionSystematics
   */
float CScission::getScissionEnergy()
{
  float R = sep + R1 + R2;
  float momInertiaTot = momInertia1 + momInertia2 + mu0*pow(R,2);
  float EE = Vc0/R + kRotate/2./momInertiaTot*pow(fJ,2);
  float mass1 = mass->getFiniteRangeMass(Z1,A1);
  float mass2 = mass->getFiniteRangeMass(Z2,A2);
  return EE + mass1 + mass2;
}
//***************************************************
  /**
   * returns the asymmetric fission scission energy
   * from a two sphere approximation using the separation sep
   * previously determined from init() to sigmaFissionSystematics()
   /param iZ1 is the proton number of one of the fission fragments
   /param iA1 is the mass number of one of the fission fragments
  */
float CScission::getScissionEnergy(int iZ1, int iA1)
{
  A1 = (float)iA1;
  Z1 = (float)iZ1;
  A2 = A-A1;
  int iA2 = iA - iA1;
  Z2 = Z - Z1;
  int iZ2 = iZ - iZ1;
  R1 = r0*pow(A1,(float)(1./3.));
  R2 = r0*pow(A2,(float)(1./3.));
  float R = sep + R1 + R2;
  momInertia1 = 0.4*A1*pow(R1,2);
  momInertia2 = 0.4*A2*pow(R2,2);
  float mu = A1*A2/(A1+A2);
  float momInertiaTot = momInertia1 + momInertia2 + mu*pow(R,2);
  float EE = Z1*Z2*e2/R + kRotate/2./momInertiaTot*pow(fJ,2);

  float massLD1 = mass->getFiniteRangeMass(iZ1,iA1);
  float massLD2 = mass->getFiniteRangeMass(iZ2,iA2);

/*  float massLD1 = mass->getLiquidDropMass(iZ1,iA1);
  float massLD2 = mass->getLiquidDropMass(iZ2,iA2);*/

  //float mass1 = mass->getExpMass(iZ1,iA1);
  //float mass2 = mass->getExpMass(iZ2,iA2);

  //Epair1 = mass->getPairing(iZ1,iA1);
  //Epair2 = mass->getPairing(iZ2,iA2);

  //Eshell1 = mass1 - Epair1 - massLD1;
  //Eshell2 = mass2 - Epair2 - massLD2;

  return  massLD1 + massLD2 + EE;
}
//****************************************
  /**
   * returns the total fission kinetic energy from a two sphere approximation
   * using the separation sep determined from either 
   * init() to sigmaFissionSystematics(), whichever was called last
   /param iZ1 is the proton number of one of the fragments
   /param iZ2 is the mass number of one of the fragments
   */
float CScission::getFissionKineticEnergy(int iZ1, int iA1)
{
  float A1 = (float)iA1;
  float A2 = A-A1;
  float Z1 = (float)iZ1;
  float Z2 = Z - Z1;
  float R1 = r0*pow(A1,(float)(1./3.));
  float R2 = r0*pow(A2,(float)(1./3.));
  float R = sep + R1 + R2;
  float momInertia1 = 0.4*A1*pow(R1,2);
  float momInertia2 = 0.4*A2*pow(R2,2);
  float mu = A1*A2/(A1+A2);
  float momInertiaOrbit = mu*pow(R,2);
  float momInertiaTot = momInertia1 + momInertia2 + momInertiaOrbit;
  float fl = fJ*momInertiaOrbit/momInertiaTot;
  EkCoul =  Z1*Z2*e2/R;
  EkRot =  kRotate/2./momInertiaOrbit*pow(fl,2); 
  ekTot = EkCoul + EkRot;
  Erotate1 = pow(fJ*momInertia1/momInertiaTot,2)*kRotate/2./momInertia1;
  Erotate2 = pow(fJ*momInertia2/momInertiaTot,2)*kRotate/2./momInertia2;

  return ekTot;
}
//**************************************************************
/**
 * estimates the standard deviation of the fission mass distributions
 * from the systematics of Rusanov et al. Physics of the Atomic Nucleus 60 
 * (1997) 683 assuming a scission-point logic.
 * subsequentally, it approximates the scission configuration as two separated
 * spheres, where the separation is adjusted to reproduce the mass distribution
\param iZ0 is the proton number
\param iA0 is the mass number
\param fJ is the angular momentum
\param fUScission is the thermal excitation energy at the scission-point in MeV
*/
float CScission::sigmaFissionSystematicsScission(int iZ0, int iA0, float fJ, 
   float fUScission)
{
  if (fUScission < 0. || fUScission > 2000)
    {
    cout << "fUScission= " << fUScission << " sigmaFissionSystematics" << endl;
    abort();
    }
  iZ = iZ0;
  iA = iA0;
  A = (float)iA;
  Z = (float)iZ;
  float A3 = pow(A,(float)(1./3.));

  // on page 684, the temp is determined with a level-density parameter 0.093A
  //we must do the same to be consistent
  float temp = sqrt(fUScission/0.093/A);

  float Z2A = pow(Z,2)/A;
  //find stiffness from Fig8c
  float d2Vdeta2;
  if (Z2A < 23.49) d2Vdeta2 = 2.105;
  else if (Z2A < 30.) d2Vdeta2 = 1.923*Z2A - 43.08;
  else if (Z2A < 33.9) d2Vdeta2 = 3.643*Z2A - 94.224;
  else d2Vdeta2 = -1.3144*Z2A + 73.45;


 //use equation 1 to get the variance from the stiffness and temp
 float sigma2 = pow(A,2)*temp/16./d2Vdeta2;   


 //correction for angular momentum from eq 17 and 18.
 float d2sdl2;
if (Z2A >32.7) d2sdl2 = -0.1310*temp - 0.05147*Z2A +0.000766*pow(Z2A,2) 
      +0.00289*temp*Z2A + .970; 
 else if (Z2A > 31.) d2sdl2 = .2873*temp + 0.03687*Z2A 
       -0.00974*temp*Z2A -1.1143;
 else d2sdl2 = 0.0111*Z2A - .334;
 if(d2sdl2<0.)
   d2sdl2 = 0.;
float correction = d2sdl2*pow(fJ,2)/2.;


//I find that use of the full correction - overestimartes the width
// so I have scaled it
// correction*= .75;

sigma2 += correction;

//now we determine d2VdA2
 float d2VdA2 = temp/sigma2;


 //this has a component of Coloumb energy of each fragment
 float alpha = Z/A;
 float Ec0 = 0.7053*pow(alpha,2)*pow(2.,1./3.)*20./9./A3;

 // from the surface energy
 float Es = -8.*pow(2.,1./3.)*17.9439*
  (1.-1.7826*pow((A-2.*Z)/A,2))/9./
   pow(A,(float)(4./3.));

 //find unaccounted
  d2VdA2 -= Ec0 + Es;


  float fact1 = pow(2.,5./3.)*A3*e2*r0*pow(alpha,2)/9.;
  float rr = pow(2.,2./3.)*A3*r0;
  float fact2 = 2.*e2*pow(alpha,2);
  float fact3 = (88.8945*A*A3 + 120.*pow(A,2) +67.5318*pow(A*A3,2) + 
		 46.3521*pow(A,3)*A3 - 76.32*pow(A,4))*pow(r0,4);
  float fact4 = (384.*A + 453.572*A*A3*A3+246.365*A*A*A3 + 64.8*pow(A,3))*
                pow(r0,3);
  float fact5 = (635.*A3*A3 + 609.562*A*A3 + 208.8*A*A)*pow(r0,2);
  float fact6 = (482.57*A3 + 288 *A)*r0;
  float fact7 = 144.;
  float fact8 = pow(A,5)*pow(r0,6)*(144.+544.286*A3*A3 + 1131.5*A*A3 + 
	       1411.2*A*A + 1167.49*pow(A*A3,2) + 579.465*pow(A,3)*A3 +
               158.184*pow(A,4));
  float fact9 = pow(A,4)*A3*A3*pow(r0,5)*(1088.57 + 3428.79*A3*A3 + 
	        5702.4*A*A3 + 5334.*A*A + 2941.9*pow(A*A3,2) + 
                730.08*pow(A,3)*A3);
  float fact10 = pow(A*r0,4)*A3*(3428.79 + 8640.*A3*A3 + 6720.42*A*A + 
				 1853.28*pow(A*A3,2));
  float fact11 = pow(A,4)*pow(r0,3)*(760. + 10885.7*A3*A3 + 9052.*A*A3 +
				     2822.4*A*A);
  float fact12 = pow(A,3)*A3*A3*pow(r0,2)*(5442.86+6857.57*A3*A3 + 
					   2851.2*A*A3);
  float fact13 = r0*(2743.03*pow(A,3)*A3 + 1728*pow(A,4));
  float fact14 = 576*pow(A,3);
  float fact15 = fJ*(fJ+1)*kRotate*64.;


  float s = 1.;
  int tries = 0;
  for(;;)
    {

      float y = fact1/pow(s+rr,2) - fact2/(s+rr);


      //contribution from spin
      float nom = fact3+fact4*s+fact5*s*s+fact6*pow(s,3)+fact7*pow(s,4);
      float denom = fact8+fact9*s+fact10*s*s+fact11*pow(s,3)+fact12*pow(s,4)+
         fact13*pow(s,5)+fact14*pow(s,6);
      float extra = fact15*nom/denom;
      y += extra;


      float dy = -2.*fact1/pow(s+rr,3) + fact2/pow(s+rr,2);
      //contribution from spin
      extra= fact15*(fact4 +2.*fact5*s+3.*fact6*s*s +4.*fact7*pow(s,3))/denom -
	fact15*nom/pow(denom,2)*(fact9+2.*fact10*s+3.*fact11*s*s+
	4.*fact12*pow(s,3) + 5.*fact13*pow(s,4)+6.*fact14*pow(s,5)); 
       dy += extra;

      float delta = y - d2VdA2;
      float deltaS = -delta/dy;
      if (fabs(delta) < .0001) break;
      s += deltaS;
      if(s<0.)
        s=1.E-3;
      else if(s>50.)
        s=50.;
      tries++;
      if (tries == 10 || isnan(s)) 
	{
	  //cout << "iZ= " << iZ << " iA= " << iA << " Usciss= " << fUScission 
	  //     << " J= " << fJ << endl;
          //cout << "sigma2= " << sigma2 << endl;
          s = 5.;
          break;
	}
    }
 


  sep = s;
  sep1 = sep;
  Esymmetric = getScissionEnergy();

 return sigma2;

}
//**************************************************************
/**
 * estimates the standard deviation of the fission mass distributions
 * from the systematics of Rusanov et al. Physics of the Atomic Nucleus 60 
 * (1997) 683 assuming a saddle-point logic.
 * subsequentally, it approximates the scission configuration as two separated
 * spheres, where the separation is adjusted to reproduce the mass distribution
\param iZ0 is the proton number
\param iA0 is the mass number
\param fJ is the angular momentum
\param fUScission is the thermal excitation energy at the scission-point in MeV
*/
float CScission::sigmaFissionSystematicsSaddle(int iZ0, int iA0, float fJ, 
   float fUScission)
{
  if (fUScission < 0. || fUScission > 2000)
    {
    cout << "fUScission= " << fUScission << " sigmaFissionSystematics" << endl;
    abort();
    }
  iZ = iZ0;
  iA = iA0;
  A = (float)iA;
  Z = (float)iZ;
  float A3 = pow(A,(float)(1./3.));

  // on page 684, the temp is determined with a level-density parameter 0.093A
  //we must do the same to be consistent
  float temp = sqrt(fUScission/0.093/A);

  float Z2A = pow(Z,2)/A;
  //find stiffness from Fig8c
  float d2Vdeta2;
  if (Z2A < 23.49) d2Vdeta2 = 2.105;
  else if (Z2A < 31.57) d2Vdeta2 = 1.923*Z2A - 43.08;
  else if (Z2A < 34.2) d2Vdeta2 = 3.19*Z2A - 83.06;
  else d2Vdeta2 = -1.7287*Z2A + 85.42;


 //use equation 1 to get the variance from the stiffness and temp
 float sigma2 = pow(A,2)*temp/16./d2Vdeta2;   


 //correction for angular momentum from eq 17 and 18.
 float d2sdl2;
if (Z2A >32.7) d2sdl2 = -0.1310*temp - 0.05147*Z2A +0.000766*pow(Z2A,2) 
      +0.00289*temp*Z2A + .970; 
 else if (Z2A > 31.) d2sdl2 = .2873*temp + 0.03687*Z2A 
       -0.00974*temp*Z2A -1.1143;
 else d2sdl2 = 0.0111*Z2A - .334;
 if(d2sdl2<0.)
   d2sdl2 = 0.;
float correction = d2sdl2*pow(fJ,2)/2.;


//I find that use of the full correction - overestimartes the width
// so I have scaled it
// correction*= .75;

sigma2 += correction;

//now we determine d2VdA2
 float d2VdA2 = temp/sigma2;


 //this has a component of Coloumb energy of each fragment
 float alpha = Z/A;
 float Ec0 = 0.7053*pow(alpha,2)*pow(2.,1./3.)*20./9./A3;

 // from the surface energy
 float Es = -8.*pow(2.,1./3.)*17.9439*
  (1.-1.7826*pow((A-2.*Z)/A,2))/9./
   pow(A,(float)(4./3.));

 //find unaccounted
  d2VdA2 -= Ec0 + Es;


  float fact1 = pow(2.,5./3.)*A3*e2*r0*pow(alpha,2)/9.;
  float rr = pow(2.,2./3.)*A3*r0;
  float fact2 = 2.*e2*pow(alpha,2);
  float fact3 = (88.8945*A*A3 + 120.*pow(A,2) +67.5318*pow(A*A3,2) + 
		 46.3521*pow(A,3)*A3 - 76.32*pow(A,4))*pow(r0,4);
  float fact4 = (384.*A + 453.572*A*A3*A3+246.365*A*A*A3 + 64.8*pow(A,3))*
                pow(r0,3);
  float fact5 = (635.*A3*A3 + 609.562*A*A3 + 208.8*A*A)*pow(r0,2);
  float fact6 = (482.57*A3 + 288 *A)*r0;
  float fact7 = 144.;
  float fact8 = pow(A,5)*pow(r0,6)*(144.+544.286*A3*A3 + 1131.5*A*A3 + 
	       1411.2*A*A + 1167.49*pow(A*A3,2) + 579.465*pow(A,3)*A3 +
               158.184*pow(A,4));
  float fact9 = pow(A,4)*A3*A3*pow(r0,5)*(1088.57 + 3428.79*A3*A3 + 
	        5702.4*A*A3 + 5334.*A*A + 2941.9*pow(A*A3,2) + 
                730.08*pow(A,3)*A3);
  float fact10 = pow(A*r0,4)*A3*(3428.79 + 8640.*A3*A3 + 6720.42*A*A + 
				 1853.28*pow(A*A3,2));
  float fact11 = pow(A,4)*pow(r0,3)*(760. + 10885.7*A3*A3 + 9052.*A*A3 +
				     2822.4*A*A);
  float fact12 = pow(A,3)*A3*A3*pow(r0,2)*(5442.86+6857.57*A3*A3 + 
					   2851.2*A*A3);
  float fact13 = r0*(2743.03*pow(A,3)*A3 + 1728*pow(A,4));
  float fact14 = 576*pow(A,3);
  float fact15 = fJ*(fJ+1)*kRotate*64.;


  float s = 1.;
  int tries = 0;
  for(;;)
    {

      float y = fact1/pow(s+rr,2) - fact2/(s+rr);


      //contribution from spin
      float nom = fact3+fact4*s+fact5*s*s+fact6*pow(s,3)+fact7*pow(s,4);
      float denom = fact8+fact9*s+fact10*s*s+fact11*pow(s,3)+fact12*pow(s,4)+
         fact13*pow(s,5)+fact14*pow(s,6);
      float extra = fact15*nom/denom;
      y += extra;


      float dy = -2.*fact1/pow(s+rr,3) + fact2/pow(s+rr,2);
      //contribution from spin
      extra= fact15*(fact4 +2.*fact5*s+3.*fact6*s*s +4.*fact7*pow(s,3))/denom -
	fact15*nom/pow(denom,2)*(fact9+2.*fact10*s+3.*fact11*s*s+
	4.*fact12*pow(s,3) + 5.*fact13*pow(s,4)+6.*fact14*pow(s,5)); 
       dy += extra;

      float delta = y - d2VdA2;
      float deltaS = -delta/dy;
      if (fabs(delta) < .0001) break;
      s += deltaS;
      if(s<0.)
        s=1.E-3;
      else if(s>50.)
        s=50.;
      tries++;
      if (tries == 10 || isnan(s)) 
	{
	  //cout << "iZ= " << iZ << " iA= " << iA << " Usciss= " << fUScission 
	  //     << " J= " << fJ << endl;
          //cout << "sigma2= " << sigma2 << endl;
          s = 5.;
          break;
	}
    }
 


  sep = s;
  sep1 = sep;
  Esymmetric = getScissionEnergy();

 return sigma2;

}
