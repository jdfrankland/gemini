#include "CMerson.h"
#include "CPotPara.h"
#include "CWaves.h"

class CScatter : public CMerson
{
 protected:
  double scale;
  double mu;
  double V;
  double VE;
  double VEE;

  double Wvol;
  double Wsur;
  double WsurE;

  double Rsur;
  double asur;

  double a;
  double zz;
  double Rc;
  double Ri;
  double ai;
  CPotPara *RealpotPara;
  CPotPara *ImagVpotPara; // volume imag
  CPotPara *ImagSpotPara; // surface imag
  static double const e2;
  static double const kconstant;

  double Rmatch;
  double gamma;
  valarray<double> diff(double, valarray<double>);
  double Kwave2;
  double Kwave;
  double muhbar;
  int l;
 public:
  double plb; //!<  pi-lambda-bar squared in mb 
  double R;
  double Rboundary;
  double rBarrier;
  double barrier;
  double omega;
  double Ecm;

  CScatter();
  CScatter(double mu0, double V0, double VE0, double R0, double a0,double zz0, double Rc0);
  virtual ~CScatter();
  void init(double mu0, double V0, double VE0, double R0, double a0,double zz0, double Rc0);
  void init(double mu0, double V0, double VE0, double VEE0, double R0, double a0,double zz0, double Rc0);
  void init(double mu0, double V0, double VE0, double VEE0, double R0, double a0,double zz0, double Rc0, double Wsur, double WsurE, double RW, double aW);

  void init(double mu0, double V0, double VE0, double R0, double a0,double zz0, double Rc0, double Wvol, double Wsur, double Ri, double ai);
  double getBoundary(int l0);
  double getBarrier(int l0);
  double getRealPotential(double r, double E);
  double getImagPotential(double r, double E);
  double getCoulPotential(double r);
  double getWaveFunction(double energyCM, int l0);
  double getTl_OM(double energyCM, int l0);

};
