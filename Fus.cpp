#include "CFus.h"
#include <cmath>

//the following is needed in the ROOT version
//ClassImp(CFus) 

float const CFus::E=3.3;
float const CFus::H=0.65;
float const CFus::D=.03;
float const CFus::G=0.0061;

/**
 * Constructor
\param plb0 is pi-lambda-squared in mb
\param dif0 is the diffuseness in hbar
*/
CFus::CFus(float plb0, float dif0)
{
  plb = plb0;
  dif = dif0;
}
/**
 * alternative constructor
\param iZprojectile is the projectile atomic number 
\param iAprojectile is projectile mass number
\param iZtarget is the target atomic number
\param iAtarget is target mass number
\param Elab is lab. energy of projectile in MeV
\param dif0 is the diffuseness in hbar
 */
CFus::CFus(int iZprojectile, int iAprojectile, 
	   int iZtarget, int iAtarget, float Elab, float dif0)
{
  iZp = iZprojectile;
  iAp = iAprojectile;
  iZt = iZtarget;
  iAt = iAtarget;
  fElab = Elab;

    // mod-TU CMass mass;
  CMass *mass = CMass::instance();   // mod-TU
  float massProjectile= mass->getExpMass(iZprojectile,iAprojectile); // mod-TU   mass now pointer
  float massTarget= mass->getExpMass(iZtarget,iAtarget);   // mod-TU   mass now pointer
  iZcn = iZprojectile + iZtarget;
  iAcn = iAprojectile + iAtarget;
  float massCN = mass->getExpMass(iZcn,iAcn);   // mod-TU   mass now pointer


  float Q=massProjectile + massTarget - massCN;
  Ecm = (float)iAtarget/(float)(iAtarget+iAprojectile)*Elab;
  vbeam = sqrt(2.*Elab/(float)iAprojectile)*.9784;
  vcm = (float)iAprojectile*vbeam/(float)(iAtarget+iAprojectile);
  Ex = Ecm + Q;
  float Ared = (float)(iAtarget*iAprojectile)/(float)(iAtarget+iAprojectile);
   plb = 656.80/Ared/Ecm;
   dif = dif0;
}

//*********************************************
/**
 * reinitialization
\param plb0 is pi-lambda-squared in mb
\param dif0 is the diffuseness in hbar
 */
void CFus::init(float plb0, float dif0)
{
  plb = plb0;
  dif = dif0;
}
//*******************************
/**
 * returns the critical amgular momentum
\param xsection is the fusion cross section in mb
*/
float CFus:: getL0(float xsection)
{
  //first get the value for dif=0
  float l00 = sqrt(xsection/plb)-1;

  float l0 = l00;
  float dl0;
  for(;;)
    {
     int il=0;
     float xsec = 0.;
     float dxsec = 0.;
     float maxExtra = 0.;
     for (;;)
       {
         float fl = (float)il;
         float trans = 1./(1.+exp((fl-l0)/dif));
         float DtransDl0 = pow(trans,2)*exp((fl-l0)/dif)/dif;
         float extra = trans*(2.*fl+1.)*plb;
         if (maxExtra < extra) maxExtra = extra;
         xsec += extra;
         dxsec += DtransDl0*(2.*fl+1.)*plb;
         if (extra < maxExtra/1000.) break;
	 il++;
       }
     float difference = xsec - xsection;
     if (fabs(difference) < .1 ) break;
     dl0 = -difference/dxsec;
     l0 += dl0;
    }

  return l0;
}
//**************************************************************
  /**
   *returns the Bass-model potential in MeV
   \param R radius in fm
   \param AL is orbital angular momentum
   */
float CFus::F(float R,float AL)
{
float S= R - R12;
return A/R+B*AL*(AL+1.)/pow(R,2)-C/(D*exp(S/E)+G*exp(S/H));
}
//*************************************************
  /**
   *returns the derivative Bass-model potential in MeV withrespect to R
   \param R radius in fm
   \param AL is orbital angular momentum
   */
float CFus::FF(float R,float AL)
{
float S=R-R12;
float x	=-A/pow(R,2)-2.*B*AL*(AL+1.)/pow(R,3);
return x+C*(D/E*exp(S/E)+G/H*exp(S/H))/pow(D*exp(S/E)+G*exp(S/H),2);
}
//********************************************************
  /**
   *returns the 2nd derivative of the Bass-model potential 
   * in MeV withrespect to R
   \param R radius in fm
   \param AL is orbital angular momentum
   */
float  CFus::FFF(float R,float AL)
{
  float  S=R-R12;
  float  x=2.*A/pow(R,3)+6.*B*AL*(AL+1.)/pow(R,4);
  float GGG=-2.*C*pow(D/E*exp(S/E)+G/H*exp(S/H),2);
  GGG=GGG/pow(D*exp(S/E)+G*exp(S/H),3);
  float HHH=C*(D/E/E*exp(S/E)+G/H/H*exp(S/H));
  HHH=HHH/pow(D*exp(S/E)+G*exp(S/H),2);
  return x + GGG + HHH;
}
//**********************************************************
  /**
   *returns the 3nd derivative of the Bass-model potential 
   * in MeV withrespect to R
   \param R radius in fm
   \param AL is orbital angular momentum
   */
float  CFus::FFFF(float R, float AL)
{
  float S=R-R12;
  float X=-6.*A/pow(R,4)-24.*B*AL*(AL+1.)/pow(R,5);
  float Y=D*exp(S/E)+G*exp(S/H);
  float Z=D/E*exp(S/E)+G/H*exp(S/H);
  float XX=D/E/E*exp(S/E)+G/H/H*exp(S/H);
  float YY=D/pow(E,3)*exp(S/E)+G/pow(H,3)*exp(S/H);
  float ZZ=-6.*C*pow(Z,3)/pow(Y,4)+6.*C*Z*XX/pow(Y,3)-C*YY/pow(Y,2);
  return X-ZZ;
}
//*******************************************************
  /**
   * returns the critical orbital angular momentum in the Bass model
   * [Nucl Phys A231 (1974) 45 ] USING THE 1977
   * BASS NUCLEAR POTENTIAL [Phys Rev Letts 39 (1977) 265 ]
   *
   */
float CFus::getBassL()
{
 const float AF=0.75;
 float AP3 = pow((float)iAp,(float)(1./3.));
 float AT3 = pow((float)iAt,(float)(1./3.));
 float AP = (float)iAp;
 float AT = (float)iAt;
 A = 1.4399*(float)(iZp*iZt);
 B=20.90*(AP+AT)/(AP*AT);
 float AL,R,DR,DDB=0.;

    float RP=1.16*AP3-1.27/AP3;
    float RT=1.16*AT3-1.27/AT3;
    int MIN=0;
    R12=RP+RT;
    U=AP*AT/(AP+AT);

    C=RP*RT/R12;

    if (FF(R12,0.) < 0.)
      {
	CL1=0.;
        return 0.;
      }
    else
      {
	CL1=sqrt(0.25+FF(R12,0.)*pow(R12,3)/(2.*B))-0.5;
       	if (FFF(R12,CL1) > 0.)
          {
//*********** DETERMINATION OF CRITICAL L AT WHICH POCKET VANISHES***
            float R=R12;
       	    MIN= (int)CL1;
            float AJ = 0.;
            float FFB = 0.;
    	    for (int J=MIN; J<=200;J++)
	      {	 
		AJ= (float)J;
                for (;;)
                  {
                  float DR=FFF(R,AJ)/FFFF(R,AJ);
        	  R=R-DR;
        	  if(fabs(DR) <= 0.00001) break;
                  }
        	if(FF(R,AJ) <0.) break;
        	FFB=FF(R,AJ);
	      }
	    CL1=AJ-1.+FFB/(FFB-FF(R,AJ));
	    //	    cout << " Potential at contact separation " << 
	    //	      "becomes repulsive at "<< CLO << " hbar\n" <<
	    //        " pocket vanishes at " << CL1 << " habr" << endl;
	  }
//********* CALCULATION OF FUSION BARRIERS FOR EACH L VALUE**
          MAX=(int) (CL1+1.);   
	  R=R12+2.;
	  for (int J=1; J<= MAX; J++)
	    {
	      AL= (float) J-1;
              for (;;)
		{
		  DR=FF(R,AL)/FFF(R,AL);
		  R=R-DR;
		  if (fabs(DR) <= 0.00001) break;
                }
	      W[J]=F(R,AL);
 	    }           
 //**********DETERMINATION OF THE CRITICAL ANGULAR MOMENTUM L1 **
	  if (MIN > 0)
            { 
	      for ( int I=MIN; I <=MAX-1; I++)
        	{
		  float AJ= (float) I;
        	  float DD=W[I+1]-F(R12,AJ);
        	  if (DD < 0.)
                    {
		      CL1=AJ-1.+DDB/(DDB-DD);
                      break;
		    }
		  DDB=DD;
		}
	    }

	  CL2=CL1/AF;
	  E1=F(R12,CL1);
	  E2=F(R12,CL2);


         
      }
    if (Ecm < W[1]) return 0.; //below the barrier
    else if (Ecm > E2) return CL2;
    else if (Ecm > E1) return sqrt(0.25+(Ecm-F(R12,0.))*R12*R12/B)-0.5;
    else
      {
       int i;
       for (i=1;i<MAX;i++)
        {
	  if (W[i] > Ecm) break;
	}
       float DEL=W[i]-W[i-1];
       return (float)(i-2)+(Ecm-W[i-1])/DEL;
      } 

}
//*********************************************************************
/**
* returns the Bass model fusion cross section in mb
*/
float CFus::getBassXsec()
{
  float L = getBassL();
  return plb*(L*L+2*L);
}
