// A transient class to associate TrackSummary with the KalRep it corresponds to.
//
// Andrei Gaponenko, 2014

#ifndef KalmanTests_inc_TrackSummaryRecoMap_hh
#define KalmanTests_inc_TrackSummaryRecoMap_hh

#include "art/Persistency/Common/Assns.h"

#include "RecoDataProducts/inc/TrackSummary.hh"
#include "KalmanTests/inc/KalRepPtrCollection.hh"

namespace mu2e {
  typedef art::Assns<KalRepPtr,TrackSummary> TrackSummaryRecoMap;
}

#endif /* KalmanTests_inc_TrackSummaryRecoMap_hh */
