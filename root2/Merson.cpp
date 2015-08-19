#include "CMerson.h"

//**************************************************************************
// merson initialization
void CMerson::initMerson(double acc0,double step_min0, double step00,int jtest0)
{
  acc = acc0;
  step_min = step_min0;
  step0 = step00;
  jtest = jtest0;
  ok = 1;
}
//**************************************************************************
// solve system of coupled ordinary differential equations
void CMerson::solveMerson(valarray<double> *y, double xstart,
   double xend)
{
  int n = y->size();
  valarray<double> yz(n);
  valarray<double> a(n);
  valarray<double> b(n);
  valarray<double> f(n);

  // rzero is a number with a magnitude of order equal to the noise level of
  //the machine i.e. in the range of the rounding errors
  double const rzero =1.0e-23;

  //jtest = test parameter related to the steplength in the following way.
  //jtest = 0, if during the calculation we get ABS(H) less than hmin (by
  //repeated halving of the steplength), than an error message is printed,
  //ok is set equal to .FALSE. followed by return to the calling program
  //jtest = 1, checking as for jtest=0, but the calculation will continue
  //with the fixed dtep length hmin


   ok = 1;

// store internally parameters in list

   double x = xstart;
   double step = step0;

   int leave = 0;
   for (;;) 
     {
       double hsv = step;
       double cof = xend - x;
       if (abs(step) >= abs(cof)) 
	 {
	   //step length is greater than remaining interval to be integrated
           step = cof;
           //if (abs(cof/hsv) < rzero) *y = w;
           leave = 1;
	   //if leave is true, then step is equal to the maximum possible
	   //steplength within the remaining part of the domain of
	   //integration.
	 }

       yz = *y;
       double ht = step/3.;

       f = diff(x,*y);
       a = ht*f;

       *y = a + yz;

       x = x + ht;
       f = diff(x,*y);

       a = 0.5*a;

       *y = 0.5*ht*f + a + yz;

       f = diff(x,*y);

       b = 4.5*ht*f;

       *y = 0.25*b + 0.75*a + yz;

       x = x + 0.5*ht;
       f = diff(x,*y);
      
       a = 2.*ht*f + a;

       *y = 3.*a - b + yz;

       x = x + 0.5*step;
       f = diff(x,*y);

      
       int accuracy = 1;

       int k = 0;
       while (accuracy && k < n)
	 {
	   b[k] = -0.5*ht*f[k] - b[k] + 2.*a[k];
	   (*y)[k] = (*y)[k] - b[k];
	   a[k] = abs(5.*acc* (*y)[k]);
	   b[k] = abs(b[k]);
	   if (abs((*y)[k]) > rzero && b[k] > a[k]) accuracy = 0;
	   k++;
	 }

       //if accuracy is false, the required accuracy for all computed values
       //was not obtained
       if (accuracy==0)
	 {
	   //halve step length      
	   cof = 0.5*step; 
           if(abs(cof) >= step_min)
	     {
	       *y = yz;
	       x = x - step;
	       step = cof;
	       leave = 0;
	     }
	   else if (jtest == 0) 
	     {
	       cout  << "**Merson error***" << endl;
               cout << "jtest= " <<jtest <<" step_min= " << step_min 
                    << " x= " << x << endl;
	       ok = 0;
               abort();
	     }
	   else
	     {
	       //continue with constant step length equal to hmin
	       step = step_min;
	       if (hsv < 0.) step = -step;
	     }
	 }
       else
	 {
	   //required accuray obtained
	   //test if step length doubling is possible
           valarray<double> g(b-0.03125*a);
	   if (g.max() <= 0.0) step = 2.0*step;

	 }
       if (leave) break;
     }

   // calculation has finished - 

}


