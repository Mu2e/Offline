#ifndef RecoDataProducts_TriggerInfo_hh
#define RecoDataProducts_TriggerInfo_hh

#include "RecoDataProducts/inc/TriggerFlag.hh"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class CaloCluster;
  class KalSeed;
  class HelixSeed;
  class TimeCluster;

  struct TriggerInfo
  {
     
    const TriggerFlag&            triggerBits() const { return _triggerBits; }
    const art::Ptr<CaloCluster>&  caloCluster() const { return _caloCluster; }
    const art::Ptr<KalSeed>&  track() const { return _track; }
    const art::Ptr<HelixSeed>&  helix() const { return _helix; }
    const art::Ptr<TimeCluster>&  hitCluster() const { return _hitCluster; }

    TriggerFlag	       _triggerBits; 
    art::Ptr<CaloCluster>  _caloCluster; 
    art::Ptr<KalSeed> _track; // associated track
    art::Ptr<HelixSeed> _helix; // associated helix
    art::Ptr<TimeCluster> _hitCluster; // associated time cluster
  };
  typedef std::vector<mu2e::TriggerInfo> TriggerInfoCollection;
} 

#endif
