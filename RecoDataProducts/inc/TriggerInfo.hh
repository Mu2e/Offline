#ifndef RecoDataProducts_TriggerInfo_hh
#define RecoDataProducts_TriggerInfo_hh

#include "canvas/Persistency/Common/Ptr.h"
#include <string>
#include <vector>

namespace mu2e {

  class CaloCluster;
  class KalSeed;
  class HelixSeed;
  class TimeCluster;
  class CaloTrigSeed;
  class CosmicTrackSeed;
  struct TriggerInfo
  {
    //accessors
    std::vector<art::Ptr<CaloCluster>>     const&  caloClusters()     const { return _caloClusters; }
    std::vector<art::Ptr<KalSeed>>         const&  tracks()           const { return _tracks; }
    std::vector<art::Ptr<HelixSeed>>       const&  helixes()          const { return _helixes; }
    std::vector<art::Ptr<TimeCluster>>     const&  hitClusters()      const { return _hitClusters; }
    std::vector<art::Ptr<CaloTrigSeed>>    const&  caloTrigSeeds()    const { return _caloTrigSeeds; }
    std::vector<art::Ptr<CosmicTrackSeed>> const&  cosmics()          const { return _cosmics; }

    //data members
    std::vector<art::Ptr<CaloCluster>>     _caloClusters;
    std::vector<art::Ptr<KalSeed>>         _tracks;        // associated tracks
    std::vector<art::Ptr<HelixSeed>>       _helixes;       // associated helices
    std::vector<art::Ptr<TimeCluster>>     _hitClusters;   // associated timeClusters
    std::vector<art::Ptr<CaloTrigSeed>>    _caloTrigSeeds; // associated CaloTrigSeeds
    std::vector<art::Ptr<CosmicTrackSeed>> _cosmics;       // associated CosmicTrackSeeds
  };
  typedef std::vector<mu2e::TriggerInfo> TriggerInfoCollection;
}

#endif
