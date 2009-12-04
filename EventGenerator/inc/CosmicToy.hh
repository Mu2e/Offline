#ifndef CosmicToy_HH
#define CosmicToy_HH
//
// A really, really, stupid model of cosmic rays.
// The purpose is to provide an example of the interface.
//
// $Id: CosmicToy.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
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

  };

} // end namespace mu2e,

#endif


