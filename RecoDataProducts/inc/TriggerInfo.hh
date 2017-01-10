#ifndef RecoDataProducts_TriggerInfo_hh
#define RecoDataProducts_TriggerInfo_hh

#include "RecoDataProducts/inc/TriggerFlag.hh"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class CaloCluster;

  class TriggerInfo
  {
      const TriggerFlag&            triggerBits() const { return _triggerBits; }
      const art::Ptr<CaloCluster>&  caloCluster() const { return _caloCluster; }


      TriggerFlag	       _triggerBits; 
      art::Ptr<CaloCluster>    _caloCluster; 
  };

} 

#endif
