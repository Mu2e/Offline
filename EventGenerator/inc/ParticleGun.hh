#ifndef EventGenerator_ParticleGun_hh
#define EventGenerator_ParticleGun_hh
//
// Shoots a single particle gun and puts its output into a generated event.
//
// $Id: ParticleGun.hh,v 1.12 2012/02/17 18:31:40 mu2ecvs Exp $
// $Author: mu2ecvs $
// $Date: 2012/02/17 18:31:40 $
//
// Original author Rob Kutschke
// Modified  by MyeongJae Lee. See docdb-2049
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
    ParticleGun(CLHEP::HepRandomEngine& engine, art::Run const& run, SimpleConfig const& config);

    // adds generated particles to the collection
    virtual void generate(GenParticleCollection& out);

  private:
    ParticleGunImpl m_gun;
  };

} // end namespace mu2e,

#endif /* EventGenerator_ParticleGun_hh */
