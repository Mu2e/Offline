#ifndef EXTMONFNALGUN_HH
#define EXTMONFNALGUN_HH

// $Id: ExtMonFNALGun.hh,v 1.1 2012/02/03 06:00:31 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/03 06:00:31 $
//
// A generator to conveniently generate particles in the acceptance of
// the ExtMonFNAl filter.  The particles are generated at the entrance
// to the filter channel, in a cone around the axis of the first
// collimator.
// 
//
// Original author Andrei Gaponenko, 2012
// 
// The position is given in the Mu2e coordinate system.
//

#include "EventGenerator/inc/GeneratorBase.hh"
#include "EventGenerator/inc/ParticleGunImpl.hh"
#include "MCDataProducts/inc/GenParticle.hh"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

// Forward references.
namespace art{ class Run; }

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class ExtMonFNALGun: public GeneratorBase {

  public:
    ExtMonFNALGun( art::Run const& run, const SimpleConfig& config );

    // adds generated particles to the collection
    virtual void generate(GenParticleCollection& out);

  private:
    ParticleGunImpl m_gun;

    CLHEP::HepRotation m_rotation;
    CLHEP::Hep3Vector  m_translation;

    // From local coordinates used by our m_gun to mu2e frame
    GenParticle transform(const GenParticle& in) const;
  };

} // end namespace mu2e,

#endif /* EXTMONFNALGUN_HH */
