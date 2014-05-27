//
// Build a dictionary.
//
// $Id: classes.h,v 1.8 2014/05/27 20:09:50 gandr Exp $
// $Author: gandr $
// $Date: 2014/05/27 20:09:50 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "art/Persistency/Common/Wrapper.h"
#include "CLHEP/Vector/ThreeVector.h"

using namespace CLHEP;
#include "KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/KalRepPtrCollection.hh"
#include "KalmanTests/inc/KalFitMC.hh"
#include "KalmanTests/inc/TrackSummaryRecoMap.hh"

template class art::Ptr<KalRep>;
template class art::Wrapper<mu2e::KalRepCollection>;
template class art::Wrapper<mu2e::KalRepPtrCollection>;
template class std::vector<mu2e::TrkStrawHitInfo>;

template class std::pair<art::Ptr<art::Ptr<KalRep> >, art::Ptr<mu2e::TrackSummary> >;
template class std::pair<art::Ptr<mu2e::TrackSummary>, art::Ptr<art::Ptr<KalRep> > >;
template class art::Assns<art::Ptr<KalRep>, mu2e::TrackSummary>;
template class art::Assns<mu2e::TrackSummary, art::Ptr<KalRep> >;
template class art::Wrapper<art::Assns<art::Ptr<KalRep>, mu2e::TrackSummary> >;
template class art::Wrapper<art::Assns<mu2e::TrackSummary, art::Ptr<KalRep> > >;
