#ifndef __TrkPatRec_PhiClusterFinder_types_hh__
#define __TrkPatRec_PhiClusterFinder_types_hh__

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

namespace art {
  class Event;
}

#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

namespace mu2e {

  namespace PhiClusterFinderTypes {

    enum { kMaxTimePeaks = 100 };

    struct Config {
      fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level")};
      fhicl::Atom<std::string> tool_type{fhicl::Name("tool_type"), fhicl::Comment("tool type: Phi Cluster Finder Diag")};
      fhicl::Atom<int> mcTruth{fhicl::Name("mcTruth"), fhicl::Comment("MC truth")};
    };

    struct Data_t {

      const art::Event*               _event;
      const TimeCluster*              _timeCluster;
      const TimeClusterCollection*    _tccol;
      TimeClusterCollection*          _tccolnew;
    };
  }

}
#endif
