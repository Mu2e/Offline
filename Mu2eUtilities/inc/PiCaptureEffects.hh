#ifndef Mu2eUtilities_PiCaptureEffects_hh
#define Mu2eUtilities_PiCaptureEffects_hh

//
//
// $Id: PiCaptureEffects.hh,v 1.1 2012/04/18 22:57:19 onoratog Exp $
// $Author: onoratog $
// $Date: 2012/04/18 22:57:19 $
//
// Original author: Gianni Onorato
//

#include<vector> 

// Mu2e includes
#include "MCDataProducts/inc/GenParticle.hh"


#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGeneral.h"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

namespace mu2e {

  class PiCaptureEffects {

  public:

    PiCaptureEffects(double probInternalConversion, double elow, double ehi, int nbins);

    ~PiCaptureEffects();

    void defineOutput();

    double internalFraction();

    bool doPhoton();

    bool doElectron();

    bool doPositron();

    GenParticle outputGamma(CLHEP::Hep3Vector pos, double time);  

    GenParticle outputElec(CLHEP::Hep3Vector pos, double time);

    GenParticle outputPosit(CLHEP::Hep3Vector pos, double time);
    
  private:

    double _probInternalConversion;
    double _elow, _ehi;
    int _nbins;
    RandomUnitSphere _randomUnitSphere;
    CLHEP::RandFlat _randFlat;
    CLHEP::RandGeneral _spectrum;
    CLHEP::RandGeneral _internalFractionalSpectrum;

    double _electMass;

    bool _doPhoton;
    double _fract;
    double _e;
    double _eElec;
    double _ePosit;
    double _elecMom;
    double _positMom;

    // Photon energy spectrum as a continuous function.
    double energySpectrum(const double x);
    
    // Compute a binned representation of the photon energy spectrum.
    std::vector<double> binnedEnergySpectrum();

    // Photon energy spectrum as a continuous function.
    double internalFractionalSpectrum(const double x);

    // Compute a binned representation of the photon energy spectrum.
    // This runs from 0 to 1; we'll take the photon energy and multiply.
    std::vector<double> internalFractionalBinnedSpectrum();
 
  };
}

#endif
