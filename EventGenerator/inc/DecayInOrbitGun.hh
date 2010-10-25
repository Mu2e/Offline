#ifndef DECAYINORBIT_HH
#define DECAYINORBIT_HH
//
// Generate some number of DIO electrons.
//
// $Id: DecayInOrbitGun.hh,v 1.6 2010/10/25 21:12:44 onoratog Exp $
// $Author: onoratog $ 
// $Date: 2010/10/25 21:12:44 $
//
// 

// Framework includes
#include "FWCore/Framework/interface/Run.h"

// Mu2e includes
#include "EventGenerator/inc/GeneratorBase.hh"
#include "EventGenerator/inc/FoilParticleGenerator.hh"

// CLHEP includes
#include "CLHEP/Random/RandPoissonQ.h"

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

    // Class for generate particles from target

    FoilParticleGenerator fGenerator;

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
    CLHEP::RandPoissonQ _randPoissonQ;

    // Diagnostic histograms.
    TH1D* _hMultiplicity;
    TH1D* _hEElec;
    TH1D* _hEElecZ;
    TH1D* _hzPosition;
    TH1D* _hcz;
    TH1D* _hphi;
    TH1D* _ht;

  };

} // end namespace mu2e,

#endif
