#ifndef MCDataProducts_PhysicalVolumeInfoCollection_hh
#define MCDataProducts_PhysicalVolumeInfoCollection_hh

//
// Define a type for a collection of PhysicalVolumeInfo objects.
//
// $Id: PhysicalVolumeInfoCollection.hh,v 1.1 2011/05/24 17:16:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:43 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "MCDataProducts/inc/PhysicalVolumeInfo.hh"

namespace mu2e {
  typedef std::vector<PhysicalVolumeInfo> PhysicalVolumeInfoCollection;
}

#endif /* MCDataProducts_PhysicalVolumeInfoCollection_hh */
