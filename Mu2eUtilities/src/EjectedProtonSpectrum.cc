// Code factored out form EjectedProtonGun.
#include "Mu2eUtilities/inc/EjectedProtonSpectrum.hh"

#include <cmath>

namespace mu2e {
  double EjectedProtonSpectrum::getWeight(double e) {

    //taken from GMC
    //
    //   Ed Hungerford  Houston University May 17 1999
    //   Rashid Djilkibaev New York University (modified) May 18 1999
    //
    //   e - proton kinetic energy (MeV)
    //   p - proton Momentum (MeV/c)
    //
    //   Generates a proton spectrum similar to that observed in
    //   u capture in Si.  JEPT 33(1971)11 and PRL 20(1967)569

    //these numbers are in MeV!!!!
    static const double emn = 1.4; // replacing par1 from GMC
    static const double par2 = 1.3279;
    static const double par3=17844.0;
    static const double par4=.32218;
    static const double par5=100.;
    static const double par6=10.014;
    static const double par7=1050.;
    static const double par8=5.103;

    double spectrumWeight;

    if (e >= 20)
      {
        spectrumWeight=par5*exp(-(e-20.)/par6);
      }

    else if(e >= 8.0 && e <= 20.0)
      {
        spectrumWeight=par7*exp(-(e-8.)/par8);
      }
    else if (e > emn)
      {
        double xw=(1.-emn/e);
        double xu=std::pow(xw,par2);
        double xv=par3*exp(-par4*e);
        spectrumWeight=xv*xu;
      }
    else
      {
        spectrumWeight = 0.;
      }
    return spectrumWeight;
  }
}
