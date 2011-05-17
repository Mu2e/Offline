#ifndef ToyDP_StrawHitCollection_hh
#define ToyDP_StrawHitCollection_hh

//
// Define a type for a collection of StrawHit objects.
//
// $Id: StrawHitCollection.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/StrawHit.hh"

namespace mu2e {
   typedef std::vector<mu2e::StrawHit> StrawHitCollection;
}

#endif /* ToyDP_StrawHitCollection_hh */
