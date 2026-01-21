#ifndef RecoDataProducts_CaloCluster_hh
#define RecoDataProducts_CaloCluster_hh
//
// Calorimeter cluster information
//
// Note: The size is set independently of the CaloHitPtrVector to work with fast clustering algorithm
//
#include "canvas/Persistency/Common/Ptr.h"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include <vector>
#include <ostream>


namespace mu2e {

   class CaloCluster {

       public:
            CaloCluster() :
               diskID_(-1),time_(0.),timeErr_(0.0),energyDep_(0.),energyDepErr_(0.),
               cog_(CLHEP::Hep3Vector(0,0,0)),caloHitsPtrVector_(),size_(0),isSplit_(false)
            {}

            CaloCluster(int iSection, float time, float timeErr, float energy, float energyErr,
                        const CLHEP::Hep3Vector& cog, CaloHitPtrVector caloHits, unsigned size, bool isSplit) :
               diskID_(iSection),time_(time),timeErr_(timeErr),energyDep_(energy),energyDepErr_(energyErr),
               cog_(cog), caloHitsPtrVector_(caloHits),size_(size),isSplit_(isSplit)
            {}

            int                       diskID()            const{return diskID_;}
            int                       size()              const{return size_;}
            float                     time()              const{return time_;}
            float                     timeErr()           const{return timeErr_;}
            float                     energyDep()         const{return energyDep_;}
            float                     energyDepErr()      const{return energyDepErr_;}
            bool                      isSplit()           const{return isSplit_;}
            const CLHEP::Hep3Vector&  cog3Vector()        const{return cog_;}
            const CaloHitPtrVector&   caloHitsPtrVector() const{return caloHitsPtrVector_;}


         private:
            int               diskID_;
            float             time_;
            float             timeErr_;
            float             energyDep_;
            float             energyDepErr_;
            CLHEP::Hep3Vector cog_;
            CaloHitPtrVector  caloHitsPtrVector_;
            unsigned          size_;
            bool              isSplit_;
   };

   using CaloClusterCollection = std::vector<mu2e::CaloCluster>;
   using CaloClusterPtrCollection = std::vector<art::Ptr<mu2e::CaloCluster>>;
}

#endif
