//
// Shoots a single particle gun and puts its output into a generated event.
//
// Original author Rob Kutschke
// Modified by MyeongJae Lee. See docdb-2049
//

#include "EventGenerator/inc/ParticleGun.hh"

#include <iostream>

// Mu2e includes
#include "DataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {

  // Conversion energy for Al.  Should come from conditions.
  static const double pEndPoint = 104.96 * CLHEP::MeV;

  ParticleGun::ParticleGun(CLHEP::HepRandomEngine& engine, art::Run const&, const SimpleConfig& config)
    : m_gun(engine,
            config.getDouble("particleGun.mean",-1.),
            static_cast<PDGCode::type>(config.getInt("particleGun.id",  PDGCode::mu_minus)),

            config.getDouble("particleGun.pmin", pEndPoint),
            config.getDouble("particleGun.pmax", pEndPoint),
            config.getString("particleGun.momentumMode", "flat"),

            RandomUnitSphereParams(config.getDouble("particleGun.czmin", -1.0),
                                   config.getDouble("particleGun.czmax",  1.0),
                                   config.getDouble("particleGun.phimin", 0. ),
                                   config.getDouble("particleGun.phimax", CLHEP::twopi)),

            config.getDouble("particleGun.tmin", 0.),
            config.getDouble("particleGun.tmax", 0.),
            config.getHep3Vector("particleGun.point", CLHEP::Hep3Vector(0.,0.,0.)),
            config.getHep3Vector("particleGun.halfLength", CLHEP::Hep3Vector(0.,0.,0.)),
            config.getString("particleGun.sourceShape", "box"),

            config.getInt("particleGun.iterationLimit", 100),
            config.getBool("particleGun.throwOnIterationLimit", false),

            config.getBool("particleGun.useDetectorCoordinateSystem", false),

            (config.getBool("particleGun.doHistograms", false) ? "ParticleGun" : ""),
            config.getBool("particleGun.doNtuples", false),

            config.getBool("particleGun.verbose",false)
            )
  {
    std::vector<double> par;
    std::vector<double> def;
    def.clear();
    config.getVectorDouble ("particleGun.momentumParameters", par, def);
    m_gun.setParameters(par);
  }

  void ParticleGun::generate( GenParticleCollection& genParts) {
    m_gun.generate(genParts);
  }

} // namespace mu2e
