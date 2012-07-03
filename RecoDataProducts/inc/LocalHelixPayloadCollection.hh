#ifndef RecoDataProducts_LocalHelixPayloadCollection_hh
#define RecoDataProducts_LocalHelixPayloadCollection_hh
//
// Define a type for a collection of LocalHelixPayload objects.
//
// $Id: LocalHelixPayloadCollection.hh,v 1.1 2012/07/03 03:27:24 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/03 03:27:24 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "RecoDataProducts/inc/LocalHelixPayload.hh"

namespace mu2e {
   typedef std::vector<mu2e::LocalHelixPayload> LocalHelixPayloadCollection;
}

#endif /* RecoDataProducts_LocalHelixPayloadCollection_hh */
