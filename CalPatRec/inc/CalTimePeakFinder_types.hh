#ifndef __CalPatRec_CalTimePeakFinder_types_hh__
#define __CalPatRec_CalTimePeakFinder_types_hh__

#include <vector>

namespace art {
  class Event;
}

namespace fhicl {
  class ParameterSet;
}

#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"

namespace mu2e {

  class KalFitResultNew;

  namespace CalTimePeakFinderTypes {
    enum { kMaxTimePeaks = 100 };
//-----------------------------------------------------------------------------
// data structure shared by CalTimePeakFinder with its plugins
//-----------------------------------------------------------------------------
    struct Data_t {
      const art::Event*               _event;
      const TimeCluster*              _timeCluster;
      const CaloClusterCollection*    ccCollection;
      const ComboHitCollection*       chcol;
      TimeClusterCollection*          _tcColl;       // 'tcColl': time cluster collection
      int                             minNHits;
      double                          minClusterEnergy;
      int                             ntc;
      std::vector<float>              dtvec;   // distance from prediction, assume < 100 time peaks
    };
  }
}
#endif
