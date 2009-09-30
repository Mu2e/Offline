#ifndef PiCapture_HH
#define PiCapture_HH
//
//
// A really, really, stupid model of photons from pi-
// capture on the targets.
//
// $Id: PiCapture.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
// 

#include "EventGenerator/inc/GeneratorBase.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

class TH1D;

namespace edm{
  class Run;
}

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class PiCapture: public GeneratorBase{

  public:
    PiCapture( edm::Run& run, const SimpleConfig& config );
    virtual ~PiCapture();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    TH1D* _piCaptureMultiplicity;
    RandomUnitSphere _randomUnitSphere;

    double _mean;

  };

} // end namespace mu2e,

#endif


