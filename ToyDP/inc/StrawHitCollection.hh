#ifndef ToyDP_StrawHitCollection_hh
#define ToyDP_StrawHitCollection_hh

//
// Define a type for a collection of StrawHit objects.
//
// $Id: StrawHitCollection.hh,v 1.1 2010/07/01 13:34:57 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/07/01 13:34:57 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/StrawHit.hh"

namespace mu2e {
   typedef std::vector<mu2e::StrawHit> StrawHitCollection;
}

#endif
