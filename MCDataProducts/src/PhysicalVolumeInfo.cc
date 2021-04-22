#include "MCDataProducts/inc/PhysicalVolumeInfo.hh"

#include <iterator>

namespace mu2e {

  //================================================================
  std::ostream& operator<<(std::ostream& ost, const PhysicalVolumeInfo& vol) {
    ost << "( "
        << vol.name() << ", "
        << vol.copyNo() <<", "
        << vol.materialName()
        << " )";
    return ost;
  }

  //================================================================
  bool operator==(const PhysicalVolumeInfo& a, const PhysicalVolumeInfo& b) {
    return
      (a.name() == b.name()) &&
      (a.copyNo() == b.copyNo()) &&
      (a.materialName() == b.materialName());
  }

  //================================================================

}
