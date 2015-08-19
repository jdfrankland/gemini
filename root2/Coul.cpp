#include "CCoul.h"
using namespace std;

//Coulomb wave functions using continued-fraction evalution of Coulomb
//Functions and their derivatives
//Journal of computational Physics 46 (1982) 171
//The values of F, G and derivative dF, and dG are correct except
//their sign may be wrong. Either all right or all the wrong sign


complex<double> CCoul::logDerH(int l, double eta, double x)
{
  double const accur = 1e-33;
  int const MaxIterations=2000;
  double const accuh = sqrt(accur);
  complex<double> one(1.,0.);
  complex<double> two(2.,0.);

  complex<double> out(0.,1. - eta/x);

  complex<double> bb(2.*(x-eta),2.);
  complex<double> aa(-pow(eta,2)- pow((double)l,2) - (double)l,eta);
  complex<double> rk(0.,0.);
  complex<double> wi(0.,2.*eta);
  complex<double> rl(0.,1/x);

  if (abSum(bb) < accuh)
    {
      rl *= aa/(aa + rk + wi);
      out +=  rl*(bb+complex<double>(0.,2.));
      aa += 2.*(rk + wi + one);
      bb += complex<double>(0.,4.);
      rk+= complex<double>(4.,0.);
    }
  complex<double> dd = one/bb;
  complex<double> dl = aa*dd*rl;
  for (;;)
    {
      complex<double>outOld = out;
      out += dl;
      rk+= two;
      aa+= rk+wi;
      bb+=complex<double>(0.,2.);
      dd = one/(aa*dd+bb);
      dl*=bb*dd-one;
      //double error = abSum(dl)/abSum(out);
      //cout << out << endl;
      if (rk.real() > 2*MaxIterations) break;
      if (abs(out-outOld) < accur) break;
    }
  out += dl;
  return out;
}
//***********************************************
double CCoul::abSum(complex<double> x)
{
  return abs(x.real()) + abs(x.imag());
}
//***********************************************
double CCoul::logDerF(int l, double eta, double x)
{
  double const tiny = 1.e-30;
  double const eps = 1.e-30;
  double l1 = (double)l;
  double l2 = l1+1.;
  double S1 = l1/x + eta/l1;
  double S2 = l2/x + eta/l2;

  double out2 = S2;
  if (out2 == 0.)out2 = tiny;
  double C = out2;
  double D = 0.;
  for (;;)
  {
   double out1 = out2;
   l++;
   l1 = (double)l;
   l2 = l1+1.;
   S1 = S2;
   S2 = l2/x + eta/l2;
   double B = S1 + S2;
   double A = -(1.+pow(eta/l1,2));
   D = B + A*D;
   if (D == 0.) D = tiny;
   C = B + A/C;
   if (C == 0.) C = tiny;
   D = 1./D;
   double delta = C*D;
   out2 = out1*delta;
   //cout << out2 << endl;
   if (abs(delta-1.) < eps)break;
  }
  return out2;
}
//******************************************
int CCoul::init(int l, double eta, double x)
{
  //double xx = eta + sqrt(pow(eta,2)+ (double)(l*(l+1)));
  //if ( x < xx) cout << " x << xx= " << xx << endl;
  double f = logDerF(l,eta,x);
  complex<double> h = logDerH(l,eta,x);
  double p = h.real();
  double q = h.imag();

  F = 1./sqrt(pow(f-p,2)/q + q);
  dF = f*F;
  G = (f-p)*F/q;
  dG = p*G - q*F;
  return 1;
}

