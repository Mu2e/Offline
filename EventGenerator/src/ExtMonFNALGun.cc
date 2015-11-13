#include "EventGenerator/inc/ExtMonFNALGun.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "DataProducts/inc/PDGCode.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "cetlib/exception.h"

#include <cmath>
#include <algorithm>

namespace mu2e {

  //================================================================
  ExtMonFNALGun::ExtMonFNALGun(art::Run const&, const SimpleConfig& config)
    : GeneratorBase()
    , m_gun(
            config.getDouble("extMonFNALGun.multiplicity",-1.),
            static_cast<PDGCode::type>(config.getInt("extMonFNALGun.pdgId")),

            config.getDouble("extMonFNALGun.pmin", GeomHandle<ExtMonFNALBuilding>()->filterMagnet().nominalMomentum()),
            config.getDouble("extMonFNALGun.pmax", GeomHandle<ExtMonFNALBuilding>()->filterMagnet().nominalMomentum()),

            RandomUnitSphereParams(-1., -cos(config.getDouble("extMonFNALGun.coneAngle")),
                                   0., 2*M_PI),

            config.getDouble("extMonFNALGun.tmin", 0.),
            config.getDouble("extMonFNALGun.tmax", 0.),

            config.getHep3Vector("extMonFNALGun.offset", CLHEP::Hep3Vector(0.,0.,0.)),
            config.getHep3Vector("extMonFNALGun.halfSize", CLHEP::Hep3Vector(0.,0.,0.)),

            (config.getBool("extMonFNALGun.doHistograms", true) ? "ExtMonFNALGun" : ""),

            config.getInt("extMonFNALGun.verbosityLevel",0)
            )
  {
    initGeom(config.getString("extMonFNALGun.reference"));
  }

  //================================================================
  ExtMonFNALGun::ExtMonFNALGun(const fhicl::ParameterSet& pset)
    : GeneratorBase()
    , m_gun(
            pset.get<double>("multiplicity"),
            PDGCode::type(pset.get<int>("pdgId")),

            pset.get<double>("pmin", GeomHandle<ExtMonFNALBuilding>()->filterMagnet().nominalMomentum()),
            pset.get<double>("pmax", GeomHandle<ExtMonFNALBuilding>()->filterMagnet().nominalMomentum()),

            RandomUnitSphereParams(-1., -cos(pset.get<double>("coneAngle")), 0., 2*M_PI),

            pset.get<double>("tmin", 0.),
            pset.get<double>("tmax", 0.),

            h3v(pset.get<std::vector<double> >("offset")),
            h3v(pset.get<std::vector<double> >("halfSize")),

            pset.get<std::string>("histDirName"),

            pset.get<int>("verbosityLevel",0)
            )
  {
    initGeom(pset.get<std::string>("reference"));
  }

  //================================================================
  void ExtMonFNALGun::initGeom(const std::string& ref) {
    if(ref == "filter") {
      m_rotation = GeomHandle<ExtMonFNALBuilding>()->collimator1RotationInMu2e();
      m_translation = GeomHandle<ExtMonFNALBuilding>()->filterEntranceInMu2e();
    }
    else if(ref == "detector") {
      m_rotation = GeomHandle<ExtMonFNAL::ExtMon>()->detectorRotationInMu2e();
      m_translation = GeomHandle<ExtMonFNAL::ExtMon>()->detectorCenterInMu2e();
    }
    else {
      throw cet::exception("BADCONFIG")
        << "Unknown value extMonFNALGun.reference == \""<<ref<<"\""
        << ".  Allowed values: filter, detector."
        << "\n";
    }
  }

  //================================================================
  CLHEP::Hep3Vector ExtMonFNALGun::h3v(const std::vector<double>& v) {
    CLHEP::Hep3Vector res;
    if(!v.empty()) {
      if(v.size() != 3) {
        throw cet::exception("BADCONFIG")
          <<"Error converting to CLHEP::Hep3Vector: input size="<<v.size()
          <<", must be 3\n";
      }
      res = CLHEP::Hep3Vector(v[0], v[1], v[2]);
    }
    return res;
  }

  //================================================================
  void ExtMonFNALGun::generate( GenParticleCollection& outParts) {
    GenParticleCollection localParts;
    m_gun.generate(localParts);

    for(GenParticleCollection::const_iterator i=localParts.begin(); i!=localParts.end(); ++i) {
      outParts.push_back(transform(*i));
    }
  }

  //================================================================
  GenParticle ExtMonFNALGun::transform(const GenParticle& in) const {
    return GenParticle(in.pdgId(),
                       GenId::extMonFNALGun,
                       m_translation + m_rotation*in.position(),
                       m_rotation * in.momentum(),
                       in.time()
                       );
  }

  //================================================================
} // namespace mu2e
