#ifndef CalPatRec_CalTimePeakFinder_types_hh
#define CalPatRec_CalTimePeakFinder_types_hh

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
    struct Config {
      fhicl::Atom<int>         diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level")};
      fhicl::Atom<std::string> tool_type{fhicl::Name("tool_type"), fhicl::Comment("tool type: CalTimePeakFinder Diag")};
      fhicl::Atom<int>         mcTruth  {fhicl::Name("mcTruth"),   fhicl::Comment("MC truth")};
    };

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
