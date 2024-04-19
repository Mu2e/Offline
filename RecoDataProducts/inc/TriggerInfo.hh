#ifndef RecoDataProducts_TriggerInfo_hh
#define RecoDataProducts_TriggerInfo_hh

#include "canvas/Persistency/Common/Ptr.h"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloTrigSeed.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include <string>
#include <vector>

namespace mu2e {

  struct TriggerInfo
  {
    //accessors
    CaloClusterCollection     const&  caloClusters()     const { return _caloClusters; }
    KalSeedCollection         const&  tracks()           const { return _tracks; }
    HelixSeedCollection       const&  helixes()          const { return _helixes; }
    TimeClusterCollection     const&  hitClusters()      const { return _hitClusters; }
    CaloTrigSeedCollection    const&  caloTrigSeeds()    const { return _caloTrigSeeds; }
    CosmicTrackSeedCollection const&  cosmics()          const { return _cosmics; }

    //data members
    CaloClusterCollection     _caloClusters;
    KalSeedCollection         _tracks;        // associated tracks
    HelixSeedCollection       _helixes;       // associated helices
    TimeClusterCollection     _hitClusters;   // associated time clusters
    CaloTrigSeedCollection    _caloTrigSeeds; // associated CaloTrigSeeds
    CosmicTrackSeedCollection _cosmics;       // associated CosmicTrackSeeds
  };
  typedef std::vector<mu2e::TriggerInfo> TriggerInfoCollection;
}

#endif
