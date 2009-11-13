#ifndef CosmicToy_HH
#define CosmicToy_HH
//
// A really, really, stupid model of cosmic rays.
// The purpose is to provide an example of the interface.
//
// $Id: CosmicToy.hh,v 1.2 2009/11/13 23:29:19 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/11/13 23:29:19 $
//
// Original author Rob Kutschke
//

#include "EventGenerator/inc/GeneratorBase.hh"

class TH1D;

namespace edm{
  class Run;
}

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class CosmicToy: public GeneratorBase{

  public:
    CosmicToy( edm::Run& run, const SimpleConfig& config );
    virtual ~CosmicToy();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    // Histogram of multiplicity.
    TH1D* _cosmicMultiplicity;

    // Mean multiplicity.
    double _mean;

    // Time range ( in ns) over which to generate events.
    double _tmin;
    double _tmax;
    double _dt;

  };

} // end namespace mu2e,

#endif


