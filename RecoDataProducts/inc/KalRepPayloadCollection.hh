#ifndef RecoDataProducts_KalRepPayloadCollection_hh
#define RecoDataProducts_KalRepPayloadCollection_hh

//
// Define a type for a collection of KalRepPayload objects.
//
// $Id: KalRepPayloadCollection.hh,v 1.1 2012/07/03 03:27:24 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/03 03:27:24 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "RecoDataProducts/inc/KalRepPayload.hh"

namespace mu2e {
   typedef std::vector<mu2e::KalRepPayload> KalRepPayloadCollection;
}

#endif /* RecoDataProducts_KalRepPayloadCollection_hh */
