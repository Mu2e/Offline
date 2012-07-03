#ifndef RecoDataProducts_KalRepExtensionPayloadCollection_hh
#define RecoDataProducts_KalRepExtensionPayloadCollection_hh

//
// Define a type for a collection of KalRepExtensionPayload objects.
//
// $Id: KalRepExtensionPayloadCollection.hh,v 1.1 2012/07/03 03:27:24 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/03 03:27:24 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "RecoDataProducts/inc/KalRepExtensionPayload.hh"

namespace mu2e {
   typedef std::vector<mu2e::KalRepExtensionPayload> KalRepExtensionPayloadCollection;
}

#endif /* RecoDataProducts_KalRepExtensionPayloadCollection_hh */
