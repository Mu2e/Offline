#include "RecoDataProducts/inc/KalSeed.hh"
#include <limits>
namespace mu2e {
  std::vector<KalSegment>::const_iterator KalSeed::nearestSegment(float fltlen)  const {
    auto retval = segments().end();
    float dmin(std::numeric_limits<float>::max());
    for(auto ikseg = segments().begin(); ikseg !=segments().end(); ikseg++) {
      if(ikseg->fmin() < fltlen && ikseg->fmax() > fltlen) {
	retval = ikseg;
	break;
      }
      float dist = std::min(fabs(ikseg->fmin()-fltlen),fabs(ikseg->fmax()-fltlen));
      if(dist < dmin){
	dmin = dist;
	retval = ikseg;
      }
    }
    return retval;
  }
}

