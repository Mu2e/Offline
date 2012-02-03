#ifndef EventGenerator_ParticleGun_hh
#define EventGenerator_ParticleGun_hh
//
// Shoots a single particle gun and puts its output into a generated event.
//
// $Id: ParticleGun.hh,v 1.11 2012/02/03 05:08:06 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/03 05:08:06 $
//
// Original author Rob Kutschke
//
// The position is given in the Mu2e coordinate system.
//

#include "EventGenerator/inc/GeneratorBase.hh"
#include "EventGenerator/inc/ParticleGunImpl.hh"

// Forward references.
namespace art{ class Run; }

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class ParticleGun: public GeneratorBase{

  public:
    ParticleGun( art::Run const& run, const SimpleConfig& config );

    // adds generated particles to the collection
    virtual void generate(GenParticleCollection& out);

  private:
    ParticleGunImpl m_gun;
  };

} // end namespace mu2e,

#endif /* EventGenerator_ParticleGun_hh */
