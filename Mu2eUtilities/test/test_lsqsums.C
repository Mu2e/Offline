///////////////////////////////////////////////////////////////////////////////
// PM: testing with ROOT:
// root [0] .L build/al9-prof-e29-p087/Offline/lib/libmu2e_Mu2eUtilities.so
// root [1] .L test_lsqsums.C
// here is how the loss of the numerical accuracy starts creeping in:
//-----------------------------------------------------------------------------
// root [56] test_line(0,0,1);
// xm: 0.00000e+00 y(0): 0.00000e+00 dydx: 1.00000e+00 yerr(0): 3.01511e-01 y(xm): 0.00000e+00 yerr(xm): 3.01511e-01
// xm: 0.00000e+00 y(0): 0.00000e+00 dydx: 1.00000e+00 yerr(0): 3.01511e-01 y(xm): 0.00000e+00 yerr(xm): 3.01511e-01
// root [57] test_line(1.e10,1.e10,1);
// xm: 1.00000e+10 y(0): 0.00000e+00 dydx: 1.00000e+00 yerr(0): 1.02423e+07 y(xm): 1.00000e+10 yerr(xm): 3.01511e-01
// xm: 1.00000e+10 y(0): 0.00000e+00 dydx: 1.00000e+00 yerr(0): 9.53463e+06 y(xm): 1.00000e+10 yerr(xm): 3.01511e-01
// root [58] test_line(1.e11,1.e11,1);
// xm: 1.00000e+11 y(0):        -nan dydx:        -nan yerr(0):         inf y(xm):        -nan yerr(xm):        -nan
// xm: 1.00000e+11 y(0): 0.00000e+00 dydx: 1.00000e+00 yerr(0): 9.53463e+07 y(xm): 1.00000e+11 yerr(xm): 3.01511e-01
//-----------------------------------------------------------------------------
// root [11] test_circle(0,0,10.);
// x0r:        0.00000 y0r:        0.00000 rr:       10.00000
// x0r:        0.00000 y0r:        0.00000 rr:       10.00000
// root [12] test_circle(100,500,10.);
// x0r:      100.00000 y0r:      500.00000 rr:       10.00000
// x0r:      100.00000 y0r:      500.00000 rr:       10.00000
// root [13] test_circle(1e7,5e7,10.);
// x0r: 11108454.68712 y0r: 54323042.25367 rr:  4462881.18449
// x0r: 10000000.00000 y0r: 50000000.00000 rr:       10.00000
///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"

//-----------------------------------------------------------------------------
int test_circle(double X0 = 0., double Y0 = 0., double R=10.) {

  LsqSums4 s4;

  double xm(0), ym(0);                  // precalculate
  int    npt(10);

  for (int i=0; i<npt-5; i++) {
    double phi = 2*M_PI/npt*i;
    double x   = X0+R*cos(phi);
    double y   = Y0+R*sin(phi);

    s4.addPoint(x,y);
    xm += x;
    ym += y;
  }
  xm /= (npt-5);
  ym /= (npt-5);

  double x0r = s4.x0();
  double y0r = s4.y0();
  double rr  = s4.radius();

  printf("x0r:%15.5f y0r:%15.5f rr:%15.5f\n",x0r,y0r,rr);
//-----------------------------------------------------------------------------
// option 2: perform calculation in the COG of points
//-----------------------------------------------------------------------------
  LsqSums4 s42(xm,ym);

  for (int i=0; i<npt-5; i++) {
    double phi = 2*M_PI/npt*i;
    double x   = X0+R*cos(phi);
    double y   = Y0+R*sin(phi);
    s42.addPoint(x,y);
  }

  x0r = s42.x0();
  y0r = s42.y0();
  rr  = s42.radius();

  printf("x0r:%15.5f y0r:%15.5f rr:%15.5f\n",x0r,y0r,rr);
  return 0;
}


//-----------------------------------------------------------------------------
int test_line(double X0 = 0., double Y0 = 0., double Slope = 1.) {
  LsqSums2 s2;

  double xm(0), ym(0);                  // precalculate
  int    npt(11);
  double step(100.);

  for (int i=0; i<npt; i++) {

    double x   = X0+step*(i-(npt-1)/2.);
    double y   = Y0+step*Slope*(i-(npt-1)/2.);

    s2.addPoint(x,y);
    xm += x;
    ym += y;
  }
  xm /= npt;
  ym /= npt;

  double y0        = s2.y0   (0);
  double dydx      = s2.dydx ();
  double y0_err    = s2.y0Err();
  double y0_xm     = s2.y0   (xm);
  double y0_xm_err = s2.y0Err(xm);

  printf("xm:%12.5e y(0):%12.5e dydx:%12.5e yerr(0):%12.5e y(xm):%12.5e yerr(xm):%12.5e\n",
         xm,y0,dydx,y0_err,y0_xm,y0_xm_err);
//-----------------------------------------------------------------------------
// option 2: perform calculation in the COG of points
//-----------------------------------------------------------------------------
  LsqSums2 s22(xm,ym);

  for (int i=0; i<npt; i++) {
    double x = X0+step*(i-(npt-1)/2.);
    double y = Y0+step*Slope*(i-(npt-1)/2.); // default slope : 0.5
    s22.addPoint(x,y);
  }

  y0        = s22.y0   (0);
  dydx      = s22.dydx ();
  y0_err    = s22.y0Err();
  y0_xm     = s22.y0   (xm);
  y0_xm_err = s22.y0Err(xm);

  printf("xm:%12.5e y(0):%12.5e dydx:%12.5e yerr(0):%12.5e y(xm):%12.5e yerr(xm):%12.5e\n",
         xm,y0,dydx,y0_err,y0_xm,y0_xm_err);

  return 0;
}
