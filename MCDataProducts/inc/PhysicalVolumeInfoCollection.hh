#ifndef MCDataProducts_PhysicalVolumeInfoCollection_hh
#define MCDataProducts_PhysicalVolumeInfoCollection_hh

//
// Define a type for a collection of PhysicalVolumeInfo objects.
//
// $Id: PhysicalVolumeInfoCollection.hh,v 1.3 2013/09/27 16:03:41 gandr Exp $
// $Author: gandr $
// $Date: 2013/09/27 16:03:41 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "MCDataProducts/inc/PhysicalVolumeInfo.hh"

namespace mu2e {
  typedef std::vector<PhysicalVolumeInfo> PhysicalVolumeInfoCollection;
}

#endif /* MCDataProducts_PhysicalVolumeInfoCollection_hh */
