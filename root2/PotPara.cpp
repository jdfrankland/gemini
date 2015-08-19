#include "CPotPara.h"
#include <iostream>

using namespace std;


//*********************************************************************
// CPotPara constructor
CPotPara::CPotPara( double V0, double R0, double a0)
  :V(V0),R(R0),a(a0)
{}
//********************************************************************
// initialization
void CPotPara::init(double V0, double R0, double a0)
{
  V = V0;
  R = R0;
  a = a0;
}
//*****************************************************************
// Wood Saxon form factor
double CPotPara::woodSaxon(double r)
{
  return -V/(1.+exp((r-R)/a));
}
//******************************************************************
//derivative of Wood Saxon *-4a
double CPotPara:: derWoodSaxon(double r)
{
  float fact = exp((r-R)/a);
  return -4.*fact*V/pow(1+fact,2);
}
//******************************************************************
//alternative surface potential - Gaussian
double CPotPara::gauss(double r)
{
  return -V*exp(-pow((r-R)/a,2));
}
//*****************************************************************
CPotPara CPotPara::copy()
{
  CPotPara out(V,R,a);
  return out;
}
//******************************************************************
void CPotPara::print()
{
  cout << "V= " << V << " R= " << R << " a= " << a << endl;
}
