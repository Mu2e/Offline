#ifndef ToyDP_ToyHitCollection_hh
#define ToyDP_ToyHitCollection_hh

//
// Define a type for a collection of ToyHits.
//
// $Id: ToyHitCollection.hh,v 1.3 2011/05/18 02:27:19 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:19 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/ToyHit.hh"

namespace mu2e {
   typedef std::vector<mu2e::ToyHit> ToyHitCollection;
}

#endif /* ToyDP_ToyHitCollection_hh */
