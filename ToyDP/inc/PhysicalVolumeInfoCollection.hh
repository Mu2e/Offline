#ifndef ToyDP_PhysicalVolumeInfoCollection_hh
#define ToyDP_PhysicalVolumeInfoCollection_hh

//
// Define a type for a collection of PhysicalVolumeInfo objects.
//
// $Id: PhysicalVolumeInfoCollection.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/PhysicalVolumeInfo.hh"

namespace mu2e {
  typedef std::vector<PhysicalVolumeInfo> PhysicalVolumeInfoCollection;
}

#endif /* ToyDP_PhysicalVolumeInfoCollection_hh */
