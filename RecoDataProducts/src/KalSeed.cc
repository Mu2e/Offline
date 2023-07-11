#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include <limits>
namespace mu2e {

  KalSeed::LHPTPtr KalSeed::loopHelixFitTrajectory() const {
    if(loopHelixFit() && segments().size() > 0){
      // initialize the piecewise trajectory with the front segment
      LHPTPtr ptraj(new KalSeed::LHPT(segments().front().loopHelix()));
      auto iseg = segments().begin(); ++iseg;
      while(iseg != segments().end()){
        if(!iseg->timeRange().null())ptraj->append(iseg->loopHelix());
        ++iseg;
      }
      return ptraj;
    } else
      return LHPTPtr();
  }

  KalSeed::CHPTPtr KalSeed::centralHelixFitTrajectory() const {
    if(centralHelixFit() && segments().size() > 0){
      // initialize the piecewise trajectory with the front segment
      CHPTPtr ptraj(new KalSeed::CHPT(segments().front().centralHelix()));
      auto iseg = segments().begin(); ++iseg;
      while(iseg != segments().end()){
        if(!iseg->timeRange().null())ptraj->append(iseg->centralHelix());
        ++iseg;
      }
      return ptraj;
    } else
      return CHPTPtr();
  }

  KalSeed::KLPTPtr KalSeed::kinematicLineFitTrajectory() const {
    if(kinematicLineFit() && segments().size() > 0){
      // initialize the piecewise trajectory with the front segment
      KLPTPtr ptraj(new KalSeed::KLPT(segments().front().kinematicLine()));
      auto iseg = segments().begin(); ++iseg;
      while(iseg != segments().end()){
        if(!iseg->timeRange().null())ptraj->append(iseg->kinematicLine());
        ++iseg;
      }
      return ptraj;
    } else
      return KLPTPtr();
  }


  std::vector<KalSegment>::const_iterator KalSeed::nearestSeg(double time)  const {
    auto retval = segments().end();
    float dmin(std::numeric_limits<float>::max());
    for(auto ikseg = segments().begin(); ikseg !=segments().end(); ikseg++) {
      if(ikseg->tmin() < time && ikseg->tmax() > time) {
        retval = ikseg;
        break;
      }
      float dist = std::min(fabs(ikseg->tmin()-time),fabs(ikseg->tmax()-time));
      if(dist < dmin){
        dmin = dist;
        retval = ikseg;
      }
    }
    return retval;
  }

  std::vector<KalSegment>::const_iterator KalSeed::nearestSegmentFlt(float fltlen)  const {
    auto retval = segments().end();
    float dmin(std::numeric_limits<float>::max());
    for(auto ikseg = segments().begin(); ikseg !=segments().end(); ikseg++) {
      // limits are kept in LOCAL flight legnth
      float lflt = ikseg->localFlt(fltlen);
      if(ikseg->fmin() < lflt && ikseg->fmax() > lflt) {
        retval = ikseg;
        break;
      }
      float dist = std::min(fabs(ikseg->fmin()-lflt),fabs(ikseg->fmax()-lflt));
      if(dist < dmin){
        dmin = dist;
        retval = ikseg;
      }
    }
    return retval;
  }

  std::vector<KalSegment>::const_iterator KalSeed::nearestSegment(float time)  const {
    auto retval = segments().end();
    double tmin(std::numeric_limits<float>::max());
    for(auto ikseg = segments().begin(); ikseg !=segments().end(); ikseg++) {
      double dt = fabs(time - 0.5*(ikseg->tmin()+ikseg->tmax()));
      if(dt < tmin){
        tmin = dt;
        retval = ikseg;
      }
    }
    return retval;
  }

  std::vector<KalSegment>::const_iterator KalSeed::nearestSegment(const XYZVectorF& pos)  const {
    auto retval = segments().end();
    double dmin(std::numeric_limits<float>::max());
    for(auto ikseg = segments().begin(); ikseg !=segments().end(); ikseg++) {
      double dist = fabs(ikseg->position3().Z()-pos.Z());
      if(dist < dmin){
        dmin = dist;
        retval = ikseg;
      }
    }
    return retval;
  }

  double KalSeed::t0Val() const {
    if(segments().size() <= 0) throw cet::exception("RECO")<<"mu2e::KalSeed: no segments" << std::endl;
// find the segment nearest z=0.  If there's one segment, we are done
    auto iseg = segments().begin();
    if(_segments.size() > 1){
      auto pvel = iseg->state().velocity();
      double vz = pvel.Z(); // to a good approximation, B is along Z
      auto pref = iseg->position3();
      double zmin = pref.Z() + vz*(iseg->tmin()-iseg->tref());
      double zmax = pref.Z() + vz*(iseg->tmax()-iseg->tref());
      if(zmin > 0.0 || zmax < 0.0){
        double mindz = std::min(fabs(zmin), fabs(zmax));
        // find the segment closest to z=0
        auto jseg = ++iseg;
        while(jseg != segments().end()) {
          pvel = jseg->state().velocity();
          vz = pvel.Z();
          pref = jseg->position3();
          zmin = pref.Z() + vz*(jseg->tmin()-jseg->tref());
          zmax = pref.Z() + vz*(jseg->tmax()-jseg->tref());
          double dz = std::min(fabs(zmin),fabs(zmax));
          if(zmin < 0.0 && zmax > 0.0){
            iseg = jseg;
            break;
          } else if(dz < mindz){
            iseg = jseg;
            mindz = dz;
          }
          ++jseg;
        }
      }
    }
    return iseg->t0Val();
  }

  HitT0 KalSeed::t0() const {
    if(segments().size() <= 0) throw cet::exception("RECO")<<"mu2e::KalSeed: no segments" << std::endl;
// find the segment nearest z=0.  If there's one segment, we are done
    auto iseg = segments().begin();
    if(_segments.size() > 1){
      auto pvel = iseg->state().velocity();
      double vz = pvel.Z(); // to a good approximation, B is along Z
      auto pref = iseg->position3();
      double zmin = pref.Z() + vz*(iseg->tmin()-iseg->tref());
      double zmax = pref.Z() + vz*(iseg->tmax()-iseg->tref());
      if(zmin > 0.0 || zmax < 0.0){
        double mindz = std::min(fabs(zmin), fabs(zmax));
        // find the segment closest to z=0
        auto jseg = ++iseg;
        while(jseg != segments().end()) {
          pvel = jseg->state().velocity();
          vz = pvel.Z();
          pref = jseg->position3();
          zmin = pref.Z() + vz*(jseg->tmin()-jseg->tref());
          zmax = pref.Z() + vz*(jseg->tmax()-jseg->tref());
          double dz = std::min(fabs(zmin),fabs(zmax));
          if(zmin < 0.0 && zmax > 0.0){
            iseg = jseg;
            break;
          } else if(dz < mindz){
            iseg = jseg;
            mindz = dz;
          }
          ++jseg;
        }
      }
    }
    return HitT0(iseg->t0Val(),1.0); //FIXME
  }
}

