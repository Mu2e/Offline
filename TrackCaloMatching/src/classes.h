//
// Build a dictionary.
//
// $Id: classes.h,v 1.3 2012/07/23 17:52:27 brownd Exp $
// $Author: brownd $
// $Date: 2012/07/23 17:52:27 $
//
// Original author G. Pezzullo
//

#include <vector>

#include "art/Persistency/Common/Wrapper.h"

#include "TrackCaloMatching/inc/TrkToCaloExtrapol.hh"
#include "KalmanTrack/KalRep.hh"
#include "TrackCaloMatching/inc/TrkToCaloExtrapolCollection.hh"
#include "TrackCaloMatching/inc/TrackClusterLink.hh"

template class art::Ptr<const KalRep * const >;

template class art::Ptr<mu2e::TrkToCaloExtrapol>;
template class std::vector<art::Ptr<mu2e::TrkToCaloExtrapol> >;

template class art::Wrapper<mu2e::TrkToCaloExtrapolCollection>;

namespace {
  struct Instantiations {
    mu2e::TrackClusterLink tcl;
  };
}
template class art::Wrapper<mu2e::TrackClusterLink>;
