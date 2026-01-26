#include "Offline/Mu2eKinKal/inc/KKFitUtilities.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
namespace mu2e {
  namespace Mu2eKinKal {
    std::shared_ptr<KinKal::SensorLine> hitLine(ComboHit const& ch, Straw const& straw,StrawResponse const& strawresponse) {
      double sprop = 2*strawresponse.halfPropV(ch.strawId(),1000.0*ch.energyDep()); // calibration uses EDep in KeV
      // construct a kinematic line trajectory from this straw. the measurement point is at the earlier signal
      KinKal::VEC3 vp0(straw.wireEnd(ch.earlyEnd()));
      KinKal::VEC3 vp1(straw.wireEnd(ch.lateEnd()));
      return std::make_shared<KinKal::SensorLine>(vp0,vp1,ch.time(),sprop);
    }

    std::shared_ptr<KinKal::SensorLine> strawLine(Straw const& straw,double time) {
      KinKal::VEC3 vp0(straw.strawEnd(StrawEnd::cal));
      KinKal::VEC3 vp1(straw.strawEnd(StrawEnd::hv));
      return std::make_shared<KinKal::SensorLine>(vp0,vp1,time,CLHEP::c_light); // time is irrelevant: use speed of light as sprop
    }
    bool inDetector(KinKal::VEC3 const& point) {
      return point.Rho() < 900.0 && fabs(point.Z()) < 1800; // numbers should come from Tracker TODO
    }
    double LorentzAngle(KinKal::ClosestApproachData const& ptca, KinKal::VEC3 const& bdir) {
      // New version to harmonize the value of "lang" and "uwirephi" in RecoDataProducts/inc/TrkStrawHitSeed
      auto ppoca = ptca.particlePoca().Vect();
      auto spoca = ptca.sensorPoca().Vect();
      auto sdir = ptca.sensorDirection();
      auto tperp = ppoca - spoca;
      auto raddir = sdir.Cross(bdir);
      if (raddir.Dot(tperp) < 0.0) raddir *= -1.0; // sign radially outwards such "tperp" and "raddir" angle is within [-pi,pi]
      auto lang = atan2(tperp.Dot(raddir),tperp.Dot(bdir));
      auto lang_folded = lang;
      if (lang > M_PI/2.0) lang_folded = M_PI - lang; // "folded" version of "lang" between [0,pi/2]
      return lang_folded;
    }
    bool insideStraw(KinKal::ClosestApproachData const& ca,Straw const& straw,double ubuffer)  {
      // compute the position along the wire and compare to the 1/2 length
      // have to translate from CLHEP, should be native to Straw TODO
      double upos = ca.sensorDirection().Dot((ca.sensorPoca().Vect() - KinKal::VEC3(straw.origin())));
      return fabs(upos) < straw.halfLength() + ubuffer;
    }

    KinKal::TimeRange timeBounds(ComboHitCollection const& chits) {
      if(chits.size() == 0) return KinKal::TimeRange();
      double tmin = std::numeric_limits<float>::max();
      double tmax = std::numeric_limits<float>::min();
      for( auto const& hit : chits) {
        // filter out hits already flagged as bad TODO
        tmin = std::min(tmin,(double)hit.correctedTime());
        tmax = std::max(tmax,(double)hit.correctedTime());
      }
      return KinKal::TimeRange(tmin,tmax);
    }

    double zMid(ComboHitCollection const& chits) {
      if(chits.size() == 0) return 0.0;
      double zmin = std::numeric_limits<float>::max();
      double zmax = std::numeric_limits<float>::min();
      for( auto const& hit : chits) {
        // filter out hits already flagged as bad TODO
        double zpos = hit.pos().z();
        zmin = std::min(zmin,zpos);
        zmax = std::max(zmax,zpos);
      }
      return 0.5*(zmin+zmax);
    }

  }
}
