#ifndef __CalPatRec_CalTimePeakFinder_types_hh__
#define __CalPatRec_CalTimePeakFinder_types_hh__

#include <vector>

namespace art {
  class Event;
};

namespace fhicl {
  class ParameterSet;
};

#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

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
