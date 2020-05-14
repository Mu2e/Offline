#ifndef RecoDataProducts_TriggerInfo_hh
#define RecoDataProducts_TriggerInfo_hh

#include "RecoDataProducts/inc/TriggerFlag.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <string>
#include <array>

namespace mu2e {

  class CaloCluster;
  class KalSeed;
  class HelixSeed;
  class TimeCluster;
  class CaloTrigSeed;
  class CosmicTrackSeed;
  struct TriggerInfo
  {
    constexpr static size_t MaxNObj = 20;

    const TriggerFlag&             triggerBits()     const { return _triggerBits; }
    const std::string&             triggerPath()     const { return _triggerPath; }
    
    //accessors
    std::array<art::Ptr<CaloCluster>,     MaxNObj> const&  caloClusters()     const { return _caloClusters; }
    std::array<art::Ptr<KalSeed>,         MaxNObj> const&  tracks()           const { return _tracks; }
    std::array<art::Ptr<HelixSeed>,       MaxNObj> const&  helixes()          const { return _helixes; }
    std::array<art::Ptr<TimeCluster>,     MaxNObj> const&  hitClusters()      const { return _hitClusters; }
    std::array<art::Ptr<CaloTrigSeed>,    MaxNObj> const&  caloTrigSeeds()    const { return _caloTrigSeeds; }
    std::array<art::Ptr<CosmicTrackSeed>, MaxNObj> const&  cosmics()          const { return _cosmics; }

    TriggerFlag	           _triggerBits{}; 
    std::string            _triggerPath; 
    std::array<art::Ptr<CaloCluster>,     MaxNObj> _caloClusters; 
    std::array<art::Ptr<KalSeed>,         MaxNObj> _tracks; // associated track
    std::array<art::Ptr<HelixSeed>,       MaxNObj> _helixes; // associated helix
    std::array<art::Ptr<TimeCluster>,     MaxNObj> _hitClusters; // associated time cluster
    std::array<art::Ptr<CaloTrigSeed>,    MaxNObj> _caloTrigSeeds; //associated CaloTrigSeed
    std::array<art::Ptr<CosmicTrackSeed>, MaxNObj> _cosmics; // associated CosmicTrackSeed
  };
  typedef std::vector<mu2e::TriggerInfo> TriggerInfoCollection;
} 

#endif
