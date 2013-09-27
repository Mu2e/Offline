#ifndef MCDataProducts_PhysicalVolumeInfoCollection_hh
#define MCDataProducts_PhysicalVolumeInfoCollection_hh

//
// Define a type for a collection of PhysicalVolumeInfo objects.
//
// $Id: PhysicalVolumeInfoCollection.hh,v 1.2 2013/09/27 14:56:14 gandr Exp $
// $Author: gandr $
// $Date: 2013/09/27 14:56:14 $
//
// Original author Rob Kutschke
//

#include "cetlib/map_vector.h"

#include "MCDataProducts/inc/PhysicalVolumeInfo.hh"

namespace mu2e {
  typedef cet::map_vector<mu2e::PhysicalVolumeInfo> PhysicalVolumeInfoCollection;
}

#endif /* MCDataProducts_PhysicalVolumeInfoCollection_hh */
