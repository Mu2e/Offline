#include "EventGenerator/inc/ExtMonFNALGun.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"
#include "DataProducts/inc/PDGCode.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "cetlib_except/exception.h"

#include <cmath>
#include <algorithm>

namespace mu2e {

  //================================================================
  ExtMonFNALGun::ExtMonFNALGun(CLHEP::HepRandomEngine& engine, const fhicl::ParameterSet& pset)
    : m_gun(engine,
            pset.get<double>("multiplicity"),
            PDGCode::type(pset.get<int>("pdgId")),

            pset.get<double>("pmin", GeomHandle<ExtMonFNALBuilding>()->filterMagnet().nominalMomentum()),
            pset.get<double>("pmax", GeomHandle<ExtMonFNALBuilding>()->filterMagnet().nominalMomentum()),

            RandomUnitSphereParams(-cos(pset.get<double>("coneAngleMin")),
                                   -cos(pset.get<double>("coneAngleMax")),
                                   0., 2*M_PI),

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
    else if(ref == "productionTargetEntrance") {
      m_rotation = GeomHandle<ProductionTarget>()->protonBeamRotation();
      m_translation = GeomHandle<ProductionTarget>()->position()
        + m_rotation*CLHEP::Hep3Vector(0., 0., GeomHandle<ProductionTarget>()->halfLength());
    }
    else if(ref == "productionTargetCenter") {
      m_rotation = GeomHandle<ProductionTarget>()->protonBeamRotation();
      m_translation = GeomHandle<ProductionTarget>()->position();
    }
    else if(ref == "productionTargetExit") {
      m_rotation = GeomHandle<ProductionTarget>()->protonBeamRotation();
      m_translation = GeomHandle<ProductionTarget>()->position()
        - m_rotation*CLHEP::Hep3Vector(0., 0., GeomHandle<ProductionTarget>()->halfLength());
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
