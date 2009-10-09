#ifndef ToyDP_CrudeStrawHitCollection_hh
#define ToyDP_CrudeStrawHitCollection_hh

//
// Define a type for a collection of CrudeStrawHit objects.
//
// $Id: CrudeStrawHitCollection.hh,v 1.1 2009/10/09 13:31:32 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/09 13:31:32 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/CrudeStrawHit.hh"

namespace mu2e {
   typedef std::vector<mu2e::CrudeStrawHit> CrudeStrawHitCollection;
}

#endif
