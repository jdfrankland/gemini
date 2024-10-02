#include "CSigCharged.h"


CSigCharged::CSigCharged(string sName0, float Zp0, float Ap0)
{

  Zp = Zp0;
  Ap = Ap0;
  if (sName0 == "neutron")
    {
      neutron = 1;
      return;
    }
  neutron = 0;
  sName = "tl/"+sName0+".inv";
   string fullName=string(GINPUT)+sName;
  ifstream ifFile (fullName.c_str());

  if (ifFile.fail() )
    {
      cout << "file " << fullName << " not found in CSigCharged" << endl;
      abort();
    }

  ifFile >> rc0 >> rc1 >> rc2;
  ifFile >> rI0 >> rI1 >> rI2;
  ifFile >> omega0 >> omega1 >> omega2 >> omega3;
  ifFile >> a0;
  ifFile >> aa0;
  ifFile.close();
  ifFile.clear();
}
//******************************************
void CSigCharged::prepare(float Z, float A)
{
  if (neutron)
    {
     n0 = 10.386*(1.-exp((-Z/24.50)));
     n1 = 1.227+0.03444*Z;
     n2 = 2.;
     barrier = .5;
     return;
    }

  float A13 = pow(A,(float)(1./3.));
  float rc = rc0/A13 + rc1 + rc2*A13;
  barrier = Zp*Z*1.44/rc;

  float rI = rc0/A13 + rI1 + rI2*A13;
  float mu = A*Ap/(A+Ap);
  InvInertia = 2.*mu*pow(rI,2)/41.563;
  omega = omega0 + omega1*A + omega2*exp(-A/omega3);

  a = a0;
  aa = aa0;

  offset = barrier/2.;
}
//*****************************************
float CSigCharged::getInverseXsec(float energy)
{

  if (neutron) return  (n0+ n1*energy)*(1.-exp(-energy/n2));
  float ddelta = sqrt(pow(offset,2)+pow(energy-barrier,2)) - offset;
  float delta;
  if (energy > barrier) delta = energy - barrier - a*ddelta;
  else delta = energy - barrier - aa*ddelta;


  float out = InvInertia*(delta + omega*log(1.+exp(-delta/omega)));
  if (out < 0.) out = 0.;
  return out;
}
//******************************************

float CSigCharged::getBarrier()
{
  return  barrier;
}
