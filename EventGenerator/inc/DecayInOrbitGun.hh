#ifndef DECAYINORBIT_HH
#define DECAYINORBIT_HH
//
// Generate some number of DIO electrons.
//
// $Id: DecayInOrbitGun.hh,v 1.4 2010/08/19 22:02:51 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/08/19 22:02:51 $
//
// For now this is limited to:
//  - Uniform over the targets.
//  - Uniform within each target.
//  - Uniform in time during the requested time interval.
//  - All of the above need to be improved at a later date.
//  - Limits on cos(theta) and phi but uniform within the range.

// Framework includes
#include "FWCore/Framework/interface/Run.h"

// Mu2e includes
#include "EventGenerator/inc/GeneratorBase.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

// CLHEP includes
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandFlat.h"

// Forward declarations outside of mu2e namespace.
class TH1D;
namespace edm {
  class Run;
}

namespace mu2e {

  // Forward declarations
  class SimpleConfig;

  class DecayInOrbitGun: public GeneratorBase{

  public:
    DecayInOrbitGun( edm::Run& run, const SimpleConfig& config );
    virtual ~DecayInOrbitGun();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    // Start: parameters that can be configured from the config file.

    // Mean number of dio electrons to generate in each event.
    // If positive, use this as the mean of a Poisson distribution.
    // If negative, generate exactly std::abs(mean) every time.
    double _mean;

    // Limits on the energy range to be generated.
    // Number of bins in the binned representation of the energy spectrum.
    double _elow;
    double _ehi;
    int    _nbins;

    // Limits on the generated direction.
    double _czmin;
    double _czmax;
    double _phimin;
    double _phimax;

    // Limits on the generated time.
    double _tmin;
    double _tmax;

    // Histogram control.
    bool _doHistograms;

    // End: parameters that can be configured from the config file.

    // Random number generators.
    RandomUnitSphere    _randomUnitSphere;
    CLHEP::RandFlat     _randFlat;
    CLHEP::RandPoissonQ _randPoissonQ;
    CLHEP::RandGeneral  _shape;

    // Diagnostic histograms.
    TH1D* _hMultiplicity;
    TH1D* _hEElec;
    TH1D* _hEElecZ;
    TH1D* _hzPosition;
    TH1D* _hcz;
    TH1D* _hphi;
    TH1D* _ht;

    // Compute the value of the energy spectrum at given energy.
    double energySpectrum(double e);

    // Build a binned representation of the energy spectrum.
    std::vector<double> binnedEnergySpectrum();

  };

} // end namespace mu2e,

#endif
