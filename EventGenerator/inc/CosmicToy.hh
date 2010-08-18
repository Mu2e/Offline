#ifndef CosmicToy_HH
#define CosmicToy_HH
//
// A really, really, stupid model of cosmic rays.
// The purpose is to provide an example of the interface.
//
// $Id: CosmicToy.hh,v 1.3 2010/08/18 22:40:15 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/08/18 22:40:15 $
//
// Original author Rob Kutschke
//

// Mu2e includes.
#include "EventGenerator/inc/GeneratorBase.hh"

// CLHEP includes
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"

// Forward references outside of namespace mu2e.
class TH1D;
namespace edm{
  class Run;
}

namespace mu2e {

  // Forward references.
  class SimpleConfig;

  class CosmicToy: public GeneratorBase{

  public:
    CosmicToy( edm::Run& run, const SimpleConfig& config );
    virtual ~CosmicToy();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    // Start: run time configurable parameters

    // Mean multiplicity. If negative, use -_mean as a fixed number
    double _mean;

    // Control making of histograms.
    bool _doHistograms;

    // end: run time configurable parameters

    // Histogram of multiplicity.
    TH1D* _hMultiplicity;
    TH1D* _hMomentum;
    TH1D* _hAngle;

    // Random number distributions.
    CLHEP::RandFlat     _randFlat;
    CLHEP::RandPoissonQ _randPoissonQ;

    // Time range ( in ns) over which to generate events.
    double _tmin;
    double _tmax;
    double _dt;

  };

} // end namespace mu2e,

#endif


