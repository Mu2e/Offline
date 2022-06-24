#ifndef RecoDataProducts_CaloProtoCluster_hh
#define RecoDataProducts_CaloProtoCluster_hh

#include "canvas/Persistency/Common/Ptr.h"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include <vector>

namespace mu2e {


   class CaloProtoCluster {

       public:
          CaloProtoCluster() :
             _time(0.),_timeErr(0.0),_energyDep(0.),_energyDepErr(0.0),_caloHitsPtrVector(),_isSplit(false)
          {}

          CaloProtoCluster(float time, float timeErr, float energy, float energyErr,
                           CaloHitPtrVector caloHit, bool isSplit) :
             _time(time),_timeErr(timeErr),_energyDep(energy),_energyDepErr(energyErr),
             _caloHitsPtrVector(caloHit),_isSplit(isSplit)
          {}

          float                    time()              const {return _time;}
          float                    timeErr()           const {return _timeErr;}
          float                    energyDep()         const {return _energyDep;}
          float                    energyDepErr()      const {return _energyDepErr;}
          const CaloHitPtrVector&  caloHitsPtrVector() const {return _caloHitsPtrVector;}
          bool                     isSplit()           const {return _isSplit;}

        private:
          float             _time;
          float             _timeErr;
          float             _energyDep;
          float             _energyDepErr;
          CaloHitPtrVector  _caloHitsPtrVector;
          bool              _isSplit;

   };

   using CaloProtoClusterCollection = std::vector<mu2e::CaloProtoCluster>;
}

#endif
