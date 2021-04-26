#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"

#include <algorithm>
#include <iterator>

namespace mu2e {
  std::ostream& operator<<(std::ostream& os, const PhysicalVolumeInfoSingleStage& coll) {
    os<<"["<<coll.size()<<"]{";
    for(const auto& e: coll) {
      os<<e.first<<"=>"<<e.second<<", ";
    }
    os<<"}";
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const PhysicalVolumeInfoMultiCollection& coll) {
    os<<"["<<coll.size()<<"]{";
    std::copy(coll.begin(), coll.end(), std::ostream_iterator<PhysicalVolumeInfoSingleStage>(os, ", "));
    os<<"}";
    return os;
  }

}
