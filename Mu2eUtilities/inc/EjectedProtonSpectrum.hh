// Kinetic energy spectrum of protons from muon capture, factored out
// from EjectedProtonGun in a form usable with BinnedSpectrum.
// More comments in the .cc file.

#ifndef Mu2eUtilities_EjectedProtonSpectrum_hh
#define Mu2eUtilities_EjectedProtonSpectrum_hh

namespace mu2e {

  // A free function can not be used with BinnedSpectrum, so make this a class.
  class EjectedProtonSpectrum {
  public:
    static double getWeight(double ek); // Input energy in MeV
  };
}

#endif/*Mu2eUtilities_EjectedProtonSpectrum_hh*/
