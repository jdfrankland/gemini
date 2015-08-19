
// The abstract base class merson is used to solve coupled first-order
// differential equations. To use one must define the differential equation
// diff which is a pure virtual function in the merson class. Otherwise their
// are two user functions. 
//               initMerson(acc,step_min,step0,jtest): initializes everything
//               solveMersion(y,xstart,xend): solves the coupled set of 
//                                             equations
//
//  xstart: start value for the domain of integration (double)
//  xend: end value for the domian of integration (double)
//
//  y:
//      TNT array containing the dependent variables. When entering the routine
//      the function values y(x),k=1,n and when returning to the calling
//      program the computed function values y(xend),k=1,n
//      the TNT array contains the number of function values n
//
//  acc: prescribed relative accuracy ( to be obtained for all function values
//      in the array y)
//      
//  step_min: absolute value of the minimum step length wanted
//  step0: initial step length
//
//jtest = test parameter related to the steplength in the following way.
//jtest = 0, if during the calculation we get ABS(step) less than step_min (by
//repeated halving of the steplength), than an error message is printed,
//ok is set equal to 0 followed by return to the calling program
//jtest = 1, checking as for jtest=0, but the calculation will continue
//with the fixed step length step_min

//                
// this class originated from the CERN LIBRARY no D 208.

// PURPOSE = step by step integration of a system of first order differential
// equations
//
// dyk(x)/dx = fk(x,y1(x),y2(x),.....,yn(x)), k=1,n
//
// with automatic error controll using the method due to merson.

//  ok:
//     a logical variable which is set equal to 1 when entering the
//     routine. the value later is depending of jtest in the following way.
//     jtest=1, ok= 1 always
//     jtest=0, ok= 0 if the steplength becomes too small ( see
//     description for jtest)
//      

//base class merson

#include <valarray>
#include <iostream>

using namespace std;

class CMerson
{
 public:
  CMerson(){}
  void initMerson (double,double,double,int);
   
  void solveMerson(valarray<double>*,double,double);
  double acc;
  double step_min;
  double step0;
  int jtest;
  int ok;
  virtual valarray<double> diff (double, valarray<double>)=0;
};
