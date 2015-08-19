#include "CScatter.h"

double const CScatter::e2=1.44;
double const CScatter::kconstant =  .048192;

CScatter::CScatter()
{
  RealpotPara = new CPotPara();
  ImagVpotPara = new CPotPara(); //volume imag
  ImagSpotPara = new CPotPara(); //surface imag
}
//**************************************************
CScatter::CScatter(double mu0, double V0, double VE0, double R0, double a0, 
double zz0, double Rc0)
{
  init(mu0,V0,VE0,R0,a0,zz0,Rc0);
  VEE = 0.;
}
//**********************************************************
CScatter::~CScatter()
{
  delete RealpotPara;
  delete ImagVpotPara;
  delete ImagSpotPara;
}
//****************************************************************************
void CScatter::init(double mu0, double V0, double VE0, double R0, double a0, 
double zz0, double Rc0)
{
  mu = mu0;
  V = V0;
  VE = VE0;
  R = R0;
  a = a0;
  zz = zz0;
  Rc = Rc0;
  RealpotPara->init(1.,R,a);
  muhbar = kconstant*mu;
  ImagVpotPara->init(0.,R,a);
  ImagSpotPara->init(1.,R,a);
  VEE = 0.;

}
//*************************************************************************
void CScatter::init(double mu0, double V0, double VE0, double VEE0, double R0, double a0, double zz0, double Rc0,double Wsur0, double WsurE0, double RW0, double aW0 )
{
  mu = mu0;
  V = V0;
  VE = VE0;
  VEE = VEE0;
  R = R0;
  a = a0;
  zz = zz0;
  Rc = Rc0;
  Wsur = Wsur0;
  WsurE = WsurE0;
  Wvol = 0.;
  Rsur = RW0;
  asur = aW0;
  RealpotPara->init(1.,R,a);
  muhbar = kconstant*mu;
  ImagVpotPara->init(1.,R,a);
  ImagSpotPara->init(1.,Rsur,asur);

}
//****************************************************************************
void CScatter::init(double mu0, double V0, double VE0, double VEE0, double R0,
 double a0, double zz0, double Rc0)
{
  init(mu0,V0,VE0,R0,a0,zz0,Rc0);

  VEE = VEE0;
}
//*************************************
void CScatter::init(double mu0, double V0, double VE0, double R0, double a0, 
		    double zz0, double Rc0,double Wvol0,double Wsur0,
                    double Ri0, double ai0)
{
  mu = mu0;
  V = V0;
  VE = VE0;
  R = R0;
  a = a0;
  zz = zz0;
  Rc = Rc0;
  RealpotPara->init(1.,R,a);
  muhbar = kconstant*mu;
  Wsur = Wsur0;
  Wvol = Wvol0;
  Ri = Ri0;
  ai = ai0;
  ImagVpotPara->init(1.,Ri,ai);
  ImagSpotPara->init(1.,Ri,ai);
}
//*************************************************************
double CScatter::getCoulPotential(double r)
{
  if (r > Rc) return e2*zz/r;
  else return e2*zz/2./Rc*(3.-pow(r/Rc,2));
}
//*************************************************************
double CScatter::getRealPotential(double r, double E)
{
  double pot =  getCoulPotential(r) + (V-VE*E-VEE*pow(E,2))*RealpotPara->woodSaxon(r);
  if (l > 0) pot +=  (double)(l*(l+1))/pow(r,2)/muhbar;
   return pot;
}
//*************************************************************
double CScatter::getImagPotential(double r, double E)
{
  double pot =  Wvol*ImagVpotPara->woodSaxon(r)
    + (Wsur-WsurE*E)*ImagSpotPara->derWoodSaxon(r);
   return pot;
}
//***********************************************************
/**
 * returns the Coulomb + centrifical barrier for the specified 
 * orbital angular momentum
 \param l is the orbital angular momentum quantum number
*/
double CScatter::getBarrier(int l0)
{
  l = l0;


  double E = 0.;
  double Uold;
  double r;
  for (int i=0;i<2;i++)
    {
     r = R + 8.;
     Uold = 0.;
     double U;
     for (;;)
       {
         U = getRealPotential(r,E);
         if (U <= Uold) break;
         Uold = U;
         r -= 0.1;
         if (r <=0.) return -1.;
       }
     E = Uold;
     if (VE == 0.) break;
    }

  rBarrier = r+.1;  

  r = rBarrier - 0.5;
  double Ubelow = getRealPotential(r,E);
  r = rBarrier + 0.5;
  double Uabove = getRealPotential(r,E);
  omega = sqrt((Uabove+Ubelow - 2.*Uold)/.25/mu)*6.466;

  return Uold;
}
//************************************************************
/**
 * returns the locations of the boundary continion - ie. the location 
 * of the minimum in the potential, if -1 returned then there is no 
 *minimum
\param l0 is the orbital angular momentum quantum number
*/
double CScatter::getBoundary(int l0)
{
  l = l0;
  double r;
  double U;
  double E = 0.;
  //find boundary condition
  if (zz == 0. && l0 == 0)
    {
      r = 5.;
      U = getRealPotential(r,E);      
    }
  else
    {
     r = R + 8.;
     bool bar = 0;
     double Uold = 0.;

     for (;;)
       {
         U = getRealPotential(r,E);
         if (U < Uold && bar == 0) //found barrier
	   {
	     bar = 1;
             barrier = Uold;
             rBarrier = r + .1;
	   }
         else if (U > Uold && bar == 1) break; //found boundary
         Uold = U;
         r -= .1;
         if (r <= 0.) return -1.;
       }
    }
  Rboundary = r + .1;
  return Rboundary;
}

//************************************************************
double CScatter::getWaveFunction(double energyCM, int l0)
{
  scale = 0.; // set imaginary potential to zero

  Ecm = energyCM;
  l = l0;
  Kwave2 = kconstant*mu*energyCM;
  Kwave = sqrt(fabs(Kwave2));

  plb = 656.80/mu/Ecm;  // pi-lambda-bar squared in mb

  //Sommerfield parameter
  gamma  = zz*e2*Kwave/energyCM/2.;  

  double Uboundary = getRealPotential(Rboundary,energyCM);
  Rmatch = R + 5.;

  // find initial and matching wave functions
  //initialise wave function
  double rhoStop = Rmatch*Kwave; 
  int lMax = max(5,l);
  CWaves outWave(rhoStop,gamma,lMax);
  //for (int i=0;i<lMax;i++) sigma[i] = outWave.sigma[i];

  // initialize the coupled differential equation solver
  initMerson(.001,.00001,.1,0);



  valarray<double> WaveFunct(4);
  WaveFunct[0] = 1.;
  WaveFunct[1] = 0.;
  if (energyCM > Uboundary)
    {
    WaveFunct[2] = 0.;
    WaveFunct[3] = -sqrt((energyCM-Uboundary)*kconstant*mu);
    }
  else
    {
     WaveFunct[2] = sqrt((Uboundary-energyCM)*kconstant*mu);
    WaveFunct[3] = 0.;
    }


  solveMerson(&WaveFunct,Rboundary,Rmatch);


  //outWave gives derivates with respect to rho = Kwave*r
  //but at this point WaveFunctOut have derivative with respect to r
  // so make them with respect to rho=kwave*r

  WaveFunct[1] /= Kwave;
  WaveFunct[3] /= Kwave;

  //cout << WaveFunct[0] << " " << WaveFunct[1] << endl;

   // match wave functions 
   //real WaveFunct = AA*F + BB*G

   double  BB = outWave.dF[l]*WaveFunct[0] 
                     - outWave.F[l]*WaveFunct[1];
   double AA = -outWave.dG[l]*WaveFunct[0] 
                      + outWave.G[l]*WaveFunct[1];

  // imaginary part => Wavefunct  = CC*F + DD*G
   double DD = outWave.dF[l]*WaveFunct[2] 
                    - outWave.F[l]*WaveFunct[3];
   double CC = -outWave.dG[l]*WaveFunct[2] 
                    + outWave.G[l]*WaveFunct[3];

   double denominator = pow(AA+DD,2) + pow(CC-BB,2);
   double etaReal = (pow(AA,2)-pow(DD,2) + pow(CC,2) 
                    - pow(BB,2))/denominator;
   double etaImag = 2.*(AA*BB+CC*DD)/denominator;


   return 1.-pow(etaReal,2)-pow(etaImag,2);
   //return 1. -(pow(BB+CC,2)+pow(AA-DD,2))/(pow(BB-CC,2)+pow(AA+DD,2));
}

//************************************************************
double CScatter::getTl_OM(double energyCM, int l0)
{
  scale = 1.; // turn on imaginary potential
  double rStart = .05;

  Ecm = energyCM;
  l = l0;
  Kwave2 = kconstant*mu*energyCM;
  Kwave = sqrt(fabs(Kwave2));

  plb = 656.80/mu/Ecm;  // pi-lambda-bar squared in mb

  //Sommerfield parameter
  gamma  = zz*e2*Kwave/energyCM/2.;  


  Rmatch = R + 5.;

  // find initial and matching wave functions
  //initialise wave function
  double rhoStop = Rmatch*Kwave; 
  int lMax = max(5,l);


  CWaves outWave(rhoStop,gamma,lMax);
  //for (int i=0;i<lMax;i++) sigma[i] = outWave.sigma[i];

  // initialize the coupled differential equation solver
  initMerson(.001,.00001,.1,0);



 // potential at start
  double Vstart = getRealPotential(rStart,Ecm);
  double Wstart = getImagPotential(rStart,Ecm);
 //derivative of potential at start
  double dVstart = (getRealPotential(rStart+.01,Ecm) - Vstart)/0.01;
  double dWstart = (getImagPotential(rStart+.01,Ecm) - Wstart)/0.01;


  valarray<double> WaveFunct(4);

  // initialize wavefunctions
  double fact = pow(rStart,l+3)/2./(double)(2*l+3);
  WaveFunct[0] = pow(rStart,l+1) 
                  - muhbar*(energyCM-Vstart)*fact; // real part
  WaveFunct[2] =  Wstart*muhbar*fact;              //imaginary part

          // derivative of wavefunction
  fact = (double)(l+3)*pow(rStart,l+2)/2./double(2*l+3);
  WaveFunct[1] = (double)(l+1)*pow(rStart,l) 
                       - muhbar*(energyCM-Vstart)*fact; // real
  WaveFunct[3] = muhbar*Wstart*fact;          // imaginary

  fact = muhbar*pow(rStart,l+3)/2./(double)(2*l+3);
  WaveFunct[1] += dVstart*fact;
  WaveFunct[3] += dWstart*fact;




  solveMerson(&WaveFunct,rStart,Rmatch);


  //outWave gives derivates with respect to rho = Kwave*r
  //but at this point WaveFunctOut have derivative with respect to r
  // so make them with respect to rho=kwave*r

  WaveFunct[1] /= Kwave;
  WaveFunct[3] /= Kwave;

  //cout << WaveFunct[0] << " " << WaveFunct[1] << endl;

   // match wave functions 
   //real WaveFunct = AA*F + BB*G

   double  BB = outWave.dF[l]*WaveFunct[0] 
                     - outWave.F[l]*WaveFunct[1];
   double AA = -outWave.dG[l]*WaveFunct[0] 
                      + outWave.G[l]*WaveFunct[1];

  // imaginary part => Wavefunct  = CC*F + DD*G
   double DD = outWave.dF[l]*WaveFunct[2] 
                    - outWave.F[l]*WaveFunct[3];
   double CC = -outWave.dG[l]*WaveFunct[2] 
                    + outWave.G[l]*WaveFunct[3];

   double denominator = pow(AA+DD,2) + pow(CC-BB,2);
   double etaReal = (pow(AA,2)-pow(DD,2) + pow(CC,2) 
                    - pow(BB,2))/denominator;
   double etaImag = 2.*(AA*BB+CC*DD)/denominator;


   return 1.-pow(etaReal,2)-pow(etaImag,2);
   //return 1. -(pow(BB+CC,2)+pow(AA-DD,2))/(pow(BB-CC,2)+pow(AA+DD,2));
}
//**************************************************************************
valarray<double> CScatter::diff(double r , valarray<double> w)
{
 //       REAL (kind=8) :: pot_real,pot_imag   
 // this subroutine is used by the mersion subroutine 
 // w[0] = REAL(u(r)) real part of  radial wave function
 //w[1] = REAL(du/dr) derrivative of real part of radial wave function
 //w[2] = IMAG(u(r)) Imaginary part of radial wave function
 //w[3] = IMAG(du/dr) derrivative of imaginary part
 //F(1:4) are the derivaties of w(1:4)




  int n = w.size();
  valarray<double> f(n);
  f[0] = w[1]; //these are equal by definition

  double potReal = getRealPotential(r,Ecm)*muhbar;
  //if ( l > 0) potReal += (double)(l*(l+1))/pow(r,2);

  f[1] = -(Kwave2 - potReal)*w[0]; 


  if (n == 4)
    {
      f[2] = w[3]; // equal by definition
      double potImag = getImagPotential(r,Ecm)*muhbar*scale;
      f[1] -= potImag*w[2];
      f[3] = -(Kwave2 - potReal)*w[2] + potImag*w[0];
    }

  return f;
}
