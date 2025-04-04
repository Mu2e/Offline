#include "Offline/EventGenerator/inc/ExtMonFNALGunImpl.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "Offline/ProductionTargetGeom/inc/ProductionTarget.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "cetlib_except/exception.h"

#include <cmath>
#include <algorithm>

namespace mu2e {

  //================================================================
  ExtMonFNALGunImpl::ExtMonFNALGunImpl(CLHEP::HepRandomEngine& engine, const Config& conf)
    : m_gun(engine,
            conf.multiplicity(),
            PDGCode::type(conf.pdgId()),

            //pset.get<double>("pmax", GeomHandle<ExtMonFNALBuilding>()->filterMagnet().nominalMomentum()),
            getpmin(conf, *GeomHandle<ExtMonFNALBuilding>()),
            getpmax(conf, *GeomHandle<ExtMonFNALBuilding>()),

            RandomUnitSphereParams(-cos(conf.coneAngleMin()),
                                   -cos(conf.coneAngleMax()),
                                   0., 2*M_PI),

            conf.tmin(),
            conf.tmax(),

            conf.offset(),
            conf.halfSize(),

            conf.histDirName(),

            conf.verbosityLevel()
            )
  {
    initGeom(conf.reference());
  }

  //================================================================
  double ExtMonFNALGunImpl::getpmin(const Config& conf, const ExtMonFNALBuilding& emb) {
    double p = emb.filter().nominalMomentum();
    conf.pmin(p);
    return p;
  }

  //================================================================
  double ExtMonFNALGunImpl::getpmax(const Config& conf, const ExtMonFNALBuilding& emb) {
    double p = emb.filter().nominalMomentum();
    conf.pmax(p);
    return p;
  }

  //================================================================
  void ExtMonFNALGunImpl::initGeom(const std::string& ref) {
    if(ref == "filter") {
      m_rotation = GeomHandle<ExtMonFNALBuilding>()->filter().collimator1().rotationInMu2e();
      m_translation = GeomHandle<ExtMonFNALBuilding>()->filter().entranceInMu2e();
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
  void ExtMonFNALGunImpl::generate( GenParticleCollection& outParts) {
    GenParticleCollection localParts;
    m_gun.generate(localParts);

    for(GenParticleCollection::const_iterator i=localParts.begin(); i!=localParts.end(); ++i) {
      outParts.push_back(transform(*i));
    }
  }

  //================================================================
  GenParticle ExtMonFNALGunImpl::transform(const GenParticle& in) const {
    return GenParticle(in.pdgId(),
                       GenId::extMonFNALGun,
                       m_translation + m_rotation*in.position(),
                       m_rotation * in.momentum(),
                       in.time()
                       );
  }

  //================================================================
} // namespace mu2e
