#include "EventGenerator/inc/ExtMonFNALGun.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNALBuilding.hh"
#include "MCDataProducts/inc/PDGCode.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"

#include <cmath>
#include <algorithm>

namespace mu2e {

  ExtMonFNALGun::ExtMonFNALGun(art::Run const&, const SimpleConfig& config)
    : GeneratorBase()
    , m_gun(
            config.getDouble("extMonFNALGun.multiplicity",-1.),
            static_cast<PDGCode::type>(config.getInt("extMonFNALGun.pdgId")),

            config.getDouble("extMonFNALGun.pmin", GeomHandle<ExtMonFNALBuilding>()->extMonFNAL_nominalMomentum()),
            config.getDouble("extMonFNALGun.pmax", GeomHandle<ExtMonFNALBuilding>()->extMonFNAL_nominalMomentum()),

            RandomUnitSphereParams(-1., -cos(config.getDouble("extMonFNALGun.coneAngle")),
                                   0., 2*M_PI),

            config.getDouble("extMonFNALGun.tmin", 0.),
            config.getDouble("extMonFNALGun.tmax", 0.),

            config.getHep3Vector("extMonFNALGun.offset", CLHEP::Hep3Vector(0.,0.,0.)),
            config.getHep3Vector("extMonFNALGun.halfSize", CLHEP::Hep3Vector(0.,0.,0.)),

            (config.getBool("extMonFNALGun.doHistograms", true) ? "ExtMonFNALGun" : ""),

            config.getBool("extMonFNALGun.verbose",false)
            )
    , m_rotation(GeomHandle<ExtMonFNALBuilding>()->collimator1RotationInMu2e())
    , m_translation(GeomHandle<ExtMonFNALBuilding>()->filterEntranceInMu2e())
  {}

  void ExtMonFNALGun::generate( GenParticleCollection& outParts) {
    GenParticleCollection localParts;
    m_gun.generate(localParts);

    for(GenParticleCollection::const_iterator i=localParts.begin(); i!=localParts.end(); ++i) {
      outParts.push_back(transform(*i));
    }
  }

  GenParticle ExtMonFNALGun::transform(const GenParticle& in) const {
    return GenParticle(in.pdgId(),
                       GenId::extMonFNALGun,
                       m_translation + m_rotation*in.position(),
                       m_rotation * in.momentum(),
                       in.time()
                       );
  }

} // namespace mu2e
