#ifndef potpara_
#define potpara_
#include <cmath>
//**********************************************************************
//give potential as a function of r
class CPotPara
{
  public:
  CPotPara(){};
  CPotPara(double,double,double);
  void init(double,double,double);  //initialization
  CPotPara copy();
  double woodSaxon (double);
  double derWoodSaxon (double);
  double gauss (double); // alternative form for surface potential
  void print();
  double V;
  double R;
  double a;
};
#endif
