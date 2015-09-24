//
// Build a dictionary.
//
// $Id: classes.h,v 1.5 2014/08/15 15:01:55 murat Exp $
// $Author: murat $
// $Date: 2014/08/15 15:01:55 $
//
// Original author G. Pezzullo
//

#include <vector>

#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/Wrapper.h"

#include "TrackCaloMatching/inc/TrkToCaloExtrapol.hh"
// #include "KalmanTrack/KalRep.hh"
#include "TrackCaloMatching/inc/TrkToCaloExtrapolCollection.hh"
#include "TrackCaloMatching/inc/TrackClusterMatch.hh"
#include "TrackCaloMatching/inc/TrackClusterMatchOld.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

/* template class art::Ptr<mu2e::CaloCluster>; */

/* template class art::Ptr<KalRep>; */

template class std::vector <mu2e::TrackClusterMatch>;
template class art::Wrapper<mu2e::TrackClusterMatchCollection>;

template class std::vector <mu2e::TrackClusterMatchOld>;
template class art::Wrapper<mu2e::TrackClusterMatchOldCollection>;

template class art::Ptr<mu2e::TrkToCaloExtrapol>;
template class std::vector<art::Ptr<mu2e::TrkToCaloExtrapol> >;
template class art::Wrapper<mu2e::TrkToCaloExtrapolCollection>;
