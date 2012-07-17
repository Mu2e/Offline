//
// Build a dictionary.
//
// $Id: classes.h,v 1.2 2012/07/17 20:03:59 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/17 20:03:59 $
//
// Original author G. Pezzullo
//

#include <vector>

#include "art/Persistency/Common/Wrapper.h"

#include "TrackCaloMatching/inc/TrkToCaloExtrapol.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrackCaloMatching/inc/TrkToCaloExtrapolCollection.hh"
#include "TrackCaloMatching/inc/TrackClusterLink.hh"

template class art::Ptr<const TrkRecoTrk * const >;

template class art::Ptr<mu2e::TrkToCaloExtrapol>;
template class std::vector<art::Ptr<mu2e::TrkToCaloExtrapol> >;

template class art::Wrapper<mu2e::TrkToCaloExtrapolCollection>;

namespace {
  struct Instantiations {
    mu2e::TrackClusterLink tcl;
  };
}
template class art::Wrapper<mu2e::TrackClusterLink>;
