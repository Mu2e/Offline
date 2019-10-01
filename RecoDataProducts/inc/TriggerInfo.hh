#ifndef RecoDataProducts_TriggerInfo_hh
#define RecoDataProducts_TriggerInfo_hh

#include "RecoDataProducts/inc/TriggerFlag.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <string>

namespace mu2e {

  class CaloCluster;
  class KalSeed;
  class HelixSeed;
  class TimeCluster;
  class CaloTrigSeed;

  struct TriggerInfo
  {
     
    const TriggerFlag&             triggerBits()     const { return _triggerBits; }
    const std::string&             triggerPath()     const { return _triggerPath; }
    const art::Ptr<CaloCluster>&   caloCluster()     const { return _caloCluster; }
    const art::Ptr<KalSeed>&       track()           const { return _track; }
    const art::Ptr<HelixSeed>&     helix()           const { return _helix; }
    const art::Ptr<TimeCluster>&   hitCluster()      const { return _hitCluster; }
    const art::Ptr<CaloTrigSeed>&  caloTrigSeed()    const { return _caloTrigSeed; }

    TriggerFlag	           _triggerBits{}; 
    std::string            _triggerPath; 
    art::Ptr<CaloCluster>  _caloCluster; 
    art::Ptr<KalSeed>      _track; // associated track
    art::Ptr<HelixSeed>    _helix; // associated helix
    art::Ptr<TimeCluster>  _hitCluster; // associated time cluster
    art::Ptr<CaloTrigSeed> _caloTrigSeed; //associated CaloTrigSeed
  };
  typedef std::vector<mu2e::TriggerInfo> TriggerInfoCollection;
} 

#endif
