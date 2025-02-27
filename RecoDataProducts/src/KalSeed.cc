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

  unsigned KalSeed::nHits(bool active) const {
    if(!active) return _hits.size();
    unsigned nactive(0);
    for(auto const& hit : hits())
      if(hit.flag().hasAllProperties(StrawHitFlag::active))++nactive;
    return nactive;
  }

  KalSeed::InterIterCol KalSeed::intersections(SurfaceId const& surfid) const{
    KalSeed::InterIterCol retval;
    for(auto iinter = _inters.begin(); iinter != _inters.end(); ++iinter) {
      if(iinter->surfid_ == surfid) retval.push_back(iinter);
    }
    return retval;
  }

  std::vector<KalSegment>::const_iterator KalSeed::nearestSegment(double time)  const {
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


  std::vector<KalSegment>::const_iterator KalSeed::t0Segment(double& t0) const {
    if(segments().size() == 0) throw cet::exception("RECO")<<"mu2e::KalSeed: no segments" << std::endl;
// start with the middle
    t0 = timeRange().mid();
    auto iseg = nearestSegment(t0);
    auto oldseg = segments().end();
    unsigned ntest(0);
    while((!iseg->timeRange().inRange(t0)) && iseg != oldseg && ntest < segments().size()){
      oldseg = iseg;
      t0 = iseg->t0Val(_status);
      iseg = nearestSegment(t0);
      ++ntest;
    }
    return iseg;
  }

  double KalSeed::t0Val() const {
    double t0;
    t0Segment(t0);
    return t0;
  }

  // deprecated legacy functions

  HitT0 KalSeed::t0() const {
    double t0;
    t0Segment(t0);
    return HitT0(t0,-1.0);
  }
}

