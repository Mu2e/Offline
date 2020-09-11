#ifndef MCDataProducts_PhysicalVolumeInfoCollection_hh
#define MCDataProducts_PhysicalVolumeInfoCollection_hh

//
// Define a type for a collection of PhysicalVolumeInfo objects.
//
//
// Original author Rob Kutschke
//

#include <vector>

#include "MCDataProducts/inc/PhysicalVolumeInfo.hh"

namespace mu2e {
  typedef std::vector<PhysicalVolumeInfo> PhysicalVolumeInfoCollection;
}

#endif /* MCDataProducts_PhysicalVolumeInfoCollection_hh */
