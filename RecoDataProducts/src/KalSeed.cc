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
}

