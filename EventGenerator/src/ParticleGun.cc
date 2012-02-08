//
// Shoots a single particle gun and puts its output into a generated event.
//
// $Id: ParticleGun.cc,v 1.18 2012/02/08 17:45:10 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/02/08 17:45:10 $
//
// Original author Rob Kutschke
//

#include "EventGenerator/inc/ParticleGun.hh"

#include <iostream>

// Mu2e includes
#include "MCDataProducts/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {

  // Conversion energy for Al.  Should come from conditions.
  static const double pEndPoint = 104.96 * CLHEP::MeV;

  ParticleGun::ParticleGun(art::Run const&, const SimpleConfig& config)
    : GeneratorBase()
    , m_gun(
            config.getDouble("particleGun.mean",-1.),
            static_cast<PDGCode::type>(config.getInt("particleGun.id",  PDGCode::mu_minus)),

            config.getDouble("particleGun.pmin", pEndPoint),
            config.getDouble("particleGun.pmax", pEndPoint),

            RandomUnitSphereParams(config.getDouble("particleGun.czmin", -1.0),
                                   config.getDouble("particleGun.czmax",  1.0),
                                   config.getDouble("particleGun.phimin", 0. ),
                                   config.getDouble("particleGun.phimax", CLHEP::twopi)),

            config.getDouble("particleGun.tmin", 0.),
            config.getDouble("particleGun.tmax", 0.),
            config.getHep3Vector("particleGun.point", CLHEP::Hep3Vector(0.,0.,0.)),
            config.getHep3Vector("particleGun.halfLength", CLHEP::Hep3Vector(0.,0.,0.)),

            (config.getBool("particleGun.doHistograms", false) ? "ParticleGun" : ""),

            config.getBool("particleGun.verbose",false)
            )
  {}

  void ParticleGun::generate( GenParticleCollection& genParts) {
    m_gun.generate(genParts);
  }

} // namespace mu2e
