#include "CWaves.h"
#include "CCoul.h"

CWaves::CWaves(double rho0, double gamma0, int lMaximum0)
  :rho(rho0),gamma(gamma0),lMaximum(lMaximum0)
{

  F = new double [lMaximum+2];
  dF = new double [lMaximum+2];
  G = new double [lMaximum+2];
  dG = new double [lMaximum+2];
  sigma = new double [lMaximum+2];

  if (gamma == 0.) sphericalBessel();
  else coulombWave();  
}
//******************************************************
// destructor
CWaves::~CWaves()
{
  delete [] F;
  delete [] dF;
  delete [] G;
  delete [] dG;
  delete [] sigma;
}
//**************************************************************

// prepares Coulomb wave functions and the phase shift
// follows B. Buck paper
void CWaves::coulombWave()
{
  // first calculate sigma - Coulomb phase shift
  // to do this, start with a calculation of the 
  // phaseshift for l=50 from asymptotic series expansion

  int const lBegin = 50;
  double lbegin1 = (double)(lBegin+1);
  double alpha = atan(gamma/lbegin1);
  double beta = sqrt(pow(gamma,2)+ pow(lbegin1,2));

  double sigma0 = alpha*((double)lBegin+1./2.) + gamma*(log(beta)-1.)
    + (-sin(alpha)/12. + sin(3.*alpha)/360./pow(beta,2) 
       - sin(5.*alpha)/1260./pow(beta,4) + sin(7.*alpha)/1680./pow(beta,6) 
       - sin(9.*alpha)/1188./pow(beta,8))/beta;

  // now use recurrence relationship to get sigma for all lower l-waves
  for (int l=lBegin-1;l>=0;l--) 
    {
      sigma0 -= atan(gamma/(double)(l+1));
      if (l <= lMaximum) sigma[l] = sigma0;
    }

  //***********************************************
  // now calculate G0
 double s = 1.;
 double t = 0.;
 double S = 0.;
 double T = 1.- gamma/rho;

 double sums = s;
 double sumS = S;
 double sumt = t;
 double sumT = T;
 for (int i=0;i<26;i++)
   {
     double denominator = 2.*rho*(double)(i+1);
     double fact1 =  (double)(2*i+1)*gamma/denominator;
     double fact2 =  (pow(gamma,2) - (double)(i*(i+1)))/denominator;
     double snew = fact1*s - fact2*t;
     double tnew = fact1*t + fact2*s;
     double Snew = fact1*S - fact2*T - snew/rho;
     double Tnew = fact1*T + fact2*S - tnew/rho;

     s = snew; t = tnew; S = Snew; T = Tnew;

     sums += s;
     sumt += t;
     sumS += S;
     sumT += T;
   }

 if (abs(sums*sumT - sumS*sumt - 1.) > 1e-4)
   {
     // cout << " Problem generating G0" << endl; 
     //cout << "s*T - S*t = " << sums*sumT - sumS*sumt << endl;
     CCoul coul;
     coul.init(0,gamma,rho);
     G[0] = coul.G;
     dG[0] = coul.dG;
     
   }
 else
   {
 double fact = -gamma*log(2.*rho) + rho + sigma[0];
 G[0] = sums*cos(fact) - sumt*sin(fact);
 dG[0] = sumS*cos(fact) - sumT*sin(fact);
   }
 //****************************************************
 // now get G1  from Eq 17 B. Buck
   G[1] = ((gamma + 1./rho)*G[0] - dG[0])/sqrt(pow(gamma,2)+1.);

 //*****************************************************
 //now use recursion relationship Eq 15 of B. Buck to get all G's
 // also Eq 17 will give us the derivatives. 

 for (int l=1;l<=lMaximum;l++)
    {
      double ll = (double)l;
      double fact1 = (2.*ll+1.)*(gamma/(ll*(ll+1.)) + 1./rho);
      double fact2 = sqrt(pow(gamma,2)+pow(ll,2))/ll;
      double fact3 = sqrt(pow(gamma,2)+pow(ll+1.,2))/(ll+1.);
      double fact4 = gamma/(ll+1) + (ll+1)/rho;
      G[l+1] = (fact1*G[l] - fact2*G[l-1])/fact3;
      dG[l] = (fact4*G[l] - fact3*G[l+1]);
    }
 //*****************************************************
 // now for the F's using recursion relationship Eq 15 of B. Buck
 // also derivative from Eq 17

 int const lstart = 60;
 double const epsilon = 1.e-60;

 double fp = 0.; // value of F/scale at l+1
 double f = epsilon; // value of F/scale at l
 double fm; // value of f at l-1

 for (int l=lstart;l>=1;l--)
   {
     double ll = (double)l;
     double fact1 = (2.*ll+1.)*(gamma/(ll*(ll+1.)) + 1./rho);
     double fact2 = sqrt(pow(gamma,2)+pow(ll,2))/ll;
     double fact3 = sqrt(pow(gamma,2)+pow(ll+1.,2))/(ll+1.);
     double fact4 = gamma/(ll+1) + (ll+1)/rho;

     fm = ( fact1*f - fact3*fp )/fact2;
     if (l-1 <= lMaximum) F[l-1] = fm;
     if (l <= lMaximum) dF[l] = fact4*f - fact3*fp;
     fp = f; f = fm;
   } 
 dF[0] = (gamma + 1./rho)*F[0] - sqrt(pow(gamma,2)+1.)*F[1];

 //***********************************************
 // finding scaling parameter for F and dF,
 // the Wronskian (Eq 16) should be unity


 double scale = dF[0]*G[0] - F[0]*dG[0];
 for (int l=0;l<=lMaximum;l++)
   {
     dF[l] /= scale;
     F[l] /= scale;
   }
}
//****************************************************

void CWaves::sphericalBessel()
{

  // start with the regular soultion

   F[0]=sin(rho)/rho;
   F[1]=(F[0]-cos(rho))/rho;

   if (lMaximum >= 2)
     {
       double sa = F[0];
       double sb = F[1];

       double f0 = 0.0;
       double f1 = 1.0E0-100;

       double f;
       int const lStart = 60;
       for (int k=lStart;k>=0;k--)
         {
           f = (2.0*(double)k+3.0)*f1/rho-f0;
           if (k <= lMaximum) F[k]=f;
           f0 = f1;
           f1 = f;
	  }
       double cs=1.;
       if (abs(sa) > abs(sb)) cs=sa/f;
       else  cs=sb/f0;
       for (int k=0;k<=lMaximum;k++)F[k]*=cs;
     }

   dF[0]=(cos(rho)-sin(rho)/rho)/rho;

   for (int k=1;k<=lMaximum;k++)dF[k]=F[k-1]-(k+1.0)*F[k]/rho;

  //************ now for the irregular
  G[0]=-cos(rho)/rho;
  G[1]=(G[0]-sin(rho))/rho;

  double f0 = G[0];
  double f1 = G[1];

  int k = 2;
  double f;  
  for (;;)
    {
      f=(2.0*(double)k-1.0)*f1/rho-f0;
      G[k]=f;
      if (abs(f) >= 1.0E+300) break;
      f0 = f1;
      f1 = f;
      k++;
      if (k > lMaximum) break;
    }

  dG[0]=(sin(rho)+cos(rho)/rho)/rho;
  for (int i=1;i<=lMaximum;i++)
    dG[i]=G[i-1]-((double)i+1.0)*G[i]/rho;

   // up until this point F,G and dF,dF are the regular and irregular
   //   spherical bessel functions and their derivative 
   //  we now need to go the equivalent form of the 
   // the regular wavefunction F and the irregular G.
  for (int k=0;k<=lMaximum;k++)
    {

      dF[k] = dF[k]*rho + F[k];
      dG[k] = -dG[k]*rho - G[k]; 
      F[k] *= rho;
      G[k] *= -rho;
 
      sigma[k] = 0.;

     if (abs(dF[k]*G[k] - F[k]*dG[k] - 1.) > .001)
       {
         if (k == 1)
	   {
	 cout << "Wronskina error" << endl;
         cout << "l= "  << k << " rho= " << rho <<  endl;
         cout << F[k] << " " << dF[k] << " " << G[k] << " " << dG[k] << endl;
         cout << dF[k]*G[k] - F[k]*dG[k] << endl;
	   }
       }
    }

}

