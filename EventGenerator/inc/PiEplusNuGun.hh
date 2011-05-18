#ifndef EventGenerator_PiEplusNuGun_hh
#define EventGenerator_PiEplusNuGun_hh

//
// Generate a positron from pi^+ -> e^+ nu
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: PiEplusNuGun.hh,v 1.6 2011/05/18 22:01:46 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 22:01:46 $
//
// For now this is limited to:
//  - Uniform over the targets.
//  - Uniform in time during the requested interval.
//  - Limits on cos(theta) and phi but uniform within the range.
//

#include "EventGenerator/inc/GeneratorBase.hh"

namespace art {
  class Run;
}

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class PiEplusNuGun: public GeneratorBase{

  public:
    PiEplusNuGun( art::Run& run, const SimpleConfig& config );
    virtual ~PiEplusNuGun();

    virtual void generate( GenParticleCollection&  );

  private:

    // Start information from the run time configuration

    // Momentum of positron from pi e nu decay.
    double _p;

    // Limits on the generated direction.
    double _czmin;
    double _czmax;
    double _phimin;
    double _phimax;

    // Limits on the generated time.
    double _tmin;
    double _tmax;

    // End information from the run time configuration.

  };

} // end namespace mu2e,

#endif /* EventGenerator_PiEplusNuGun_hh */


