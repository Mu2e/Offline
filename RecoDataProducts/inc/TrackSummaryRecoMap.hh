// A transient class to associate TrackSummary with the KalRep it corresponds to.
//
// Andrei Gaponenko, 2014

#ifndef RecoDataProducts_inc_TrackSummaryRecoMap_hh
#define RecoDataProducts_inc_TrackSummaryRecoMap_hh

#include "canvas/Persistency/Common/Assns.h"

#include "Offline/RecoDataProducts/inc/TrackSummary.hh"
#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"

namespace mu2e {
  typedef art::Assns<KalRepPtr,TrackSummary> TrackSummaryRecoMap;
}

#endif /* RecoDataProducts_inc_TrackSummaryRecoMap_hh */
