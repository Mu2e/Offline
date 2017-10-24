#ifndef __CalPatRec_CalTimePeakFinder_types_hh__
#define __CalPatRec_CalTimePeakFinder_types_hh__

#include <vector>

class TH1F;
class TH2F;

namespace art {
  class Event;
};

namespace fhicl {
  class ParameterSet;
};

#include "RecoDataProducts/inc/TimeCluster.hh"

namespace mu2e {

  class KalFitResultNew;
  class CaloCluster;
  
#ifndef __CalPatRec_CalTimePeak_hh__
  class CalTimePeakCollection;
#endif
  
  namespace CalTimePeakFinderTypes {
//-----------------------------------------------------------------------------
// data structure shared by CalTimePeakFinder with its plugins
//-----------------------------------------------------------------------------
    struct Data_t {                        
      const art::Event*      _event;
      const TimeCluster*     _timeCluster;
      CalTimePeakCollection* _tpeaks;      // cache of time peaks
      TimeClusterCollection* _outseeds;
      int                    _minNHits;
    };
  }
}
#endif
