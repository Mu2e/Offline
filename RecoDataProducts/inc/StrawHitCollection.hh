#ifndef RecoDataProducts_StrawHitCollection_hh
#define RecoDataProducts_StrawHitCollection_hh

//
// Define a type for a collection of StrawHit objects.
//
// $Id: StrawHitCollection.hh,v 1.1 2011/05/24 17:16:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:43 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "RecoDataProducts/inc/StrawHit.hh"

namespace mu2e {
   typedef std::vector<mu2e::StrawHit> StrawHitCollection;
}

#endif /* RecoDataProducts_StrawHitCollection_hh */
