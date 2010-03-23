#ifndef ToyDP_PhysicalVolumeInfoCollection_hh
#define ToyDP_PhysicalVolumeInfoCollection_hh

//
// Define a type for a collection of PhysicalVolumeInfo objects.
//
// $Id: PhysicalVolumeInfoCollection.hh,v 1.1 2010/03/23 20:34:30 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/03/23 20:34:30 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/PhysicalVolumeInfo.hh"

namespace mu2e {
  typedef std::vector<PhysicalVolumeInfo> PhysicalVolumeInfoCollection;
}

#endif
