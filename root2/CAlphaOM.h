#include "CScatter.h"
#include "CAvrigeanu.h"

class CAlphaOM
{
 private:
  //CAlpha alpha;
  CAvrigeanu alpha;
  CScatter scatter;
  double sl[200]; //!< reaction cross section for each ell value
  double prob[200]; //!< cumulative probability for ell value
  double xsec; //!< total reaction cross section
  int lmax; //!< maximum ell wave considered
 public:

  CAlphaOM(double iZtar, double Atar, double Elab); 
  double getSigL(int l); 
  double getXsec();
  int getLmax();
  int getL(double random);

};
