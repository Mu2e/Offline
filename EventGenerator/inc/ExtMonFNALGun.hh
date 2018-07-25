#ifndef EXTMONFNALGUN_HH
#define EXTMONFNALGUN_HH

// $Id: ExtMonFNALGun.hh,v 1.2 2012/11/01 23:44:24 gandr Exp $
// $Author: gandr $
// $Date: 2012/11/01 23:44:24 $
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

#include "fhiclcpp/ParameterSet.h"

#include "EventGenerator/inc/ParticleGunImpl.hh"
#include "MCDataProducts/inc/GenParticle.hh"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

// Forward references.
namespace art{ class Run; }

namespace mu2e {

  class ExtMonFNALGun {

  public:
    explicit ExtMonFNALGun(CLHEP::HepRandomEngine& engine, const fhicl::ParameterSet& pset);

    // adds generated particles to the collection
    virtual void generate(GenParticleCollection& out);

  private:
    ParticleGunImpl m_gun;

    CLHEP::HepRotation m_rotation;
    CLHEP::Hep3Vector  m_translation;

    void initGeom(const std::string& refPoint);

    // From local coordinates used by our m_gun to mu2e frame
    GenParticle transform(const GenParticle& in) const;

    CLHEP::Hep3Vector h3v(const std::vector<double>& v);
  };

} // end namespace mu2e,

#endif /* EXTMONFNALGUN_HH */
