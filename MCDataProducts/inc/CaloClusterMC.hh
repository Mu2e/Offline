//
// This calss contains the truth matched CaloDigis for a calo cluster
//
// For convenience, energyDeposits() returns a vector of CaloEDepMC for each simParticle linked to a CaloCluster, 
// built from the CaloDigis (taking care of overalpping digis). 
// THE USER TAKES OWNERSHIP OF THIS OBJECT
//

#ifndef MCDataProducts_CaloClusterMC_hh
#define MCDataProducts_CaloClusterMC_hh

#include "MCDataProducts/inc/CaloHitMC.hh"
#include "MCDataProducts/inc/CaloEDepMC.hh"
#include <vector>

namespace mu2e
{   
   class CaloClusterMC
   {
       public:                    
          using CaloHitMCPtr = art::Ptr<CaloHitMC>;

          CaloClusterMC()                                        : digis_()      {};
          CaloClusterMC(const std::vector<CaloHitMCPtr>& digis) : digis_(digis) {};
                
          const std::vector<CaloHitMCPtr>&  caloDigiMCs     () const {return digis_; }
          std::vector<CaloEDepMC>            energyDeposits  () const;
          float                              totalEnergyDep  () const;
          float                              totalEnergyDepG4() const;

       private:
          std::vector<CaloHitMCPtr> digis_; 
   };
      
   using  CaloClusterMCCollection = std::vector<mu2e::CaloClusterMC>;
} 

#endif
