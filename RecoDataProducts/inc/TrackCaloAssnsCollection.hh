#ifndef TrackCaloAssns_TrackCaloAssnsCollection_hh
#define TrackCaloAssns_TrackCaloAssnsCollection_hh

//
// Define a type for a collection of TrackCaloAssns objects.
//
//
// Original author Vadim Rusu
//

#include <vector>

#include "RecoDataProducts/inc/TrackCaloAssns.hh"

namespace mu2e {
   typedef std::vector<mu2e::TrackCaloMatchInfo> TrackCaloAssnsCollection;
}

#endif 
