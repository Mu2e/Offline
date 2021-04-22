#include "RecoDataProducts/inc/KalSeed.hh"
#include <limits>
namespace mu2e {
  std::vector<KalSegment>::const_iterator KalSeed::nearestSegment(float fltlen)  const {
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

  std::vector<KalSegment>::const_iterator KalSeed::nearestSegment(const XYZVec& pos)  const {
    auto retval = segments().end();
    float dmin(std::numeric_limits<float>::max());
    for(auto ikseg = segments().begin(); ikseg !=segments().end(); ikseg++) {
      // limits are kept in LOCAL flight legnth
      float fltlen = 0.0;
      ikseg->helix().position(pos, fltlen); // calculate the fltlen from pos
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
    return iseg->t0();
  }
}

