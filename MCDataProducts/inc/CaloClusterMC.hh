#ifndef MCDataProducts_CaloClusterMC_hh
#define MCDataProducts_CaloClusterMC_hh

#include "MCDataProducts/inc/CaloEDepMC.hh"
#include <vector>

namespace mu2e
{   
   
   class CaloClusterMC
   {
       public:                    
          CaloClusterMC()                                     : edeps_()      {};
          CaloClusterMC(const std::vector<CaloEDepMC>& edeps) : edeps_(edeps) {};
                
                std::vector<CaloEDepMC>&   energyDeposits  ()                 {return edeps_; }
                CaloEDepMC&                energyDeposit   (unsigned i)       {return edeps_.at(i); } 
          const std::vector<CaloEDepMC>&   energyDeposits  ()           const {return edeps_; }
          const CaloEDepMC&                energyDeposit   (unsigned i) const {return edeps_.at(i); } 
                unsigned                   nParticles      ()           const {return edeps_.size();}
                float                      time            ()           const {return edeps_.at(0).time();}
                float                      totalEnergyDep  ()           const;
                float                      totalEnergyDepG4()           const;
                bool                       isConversion    ()           const;
                bool                       isOnlyConversion()           const;

       private:
          std::vector<CaloEDepMC> edeps_; 
   };
   
   
   typedef std::vector<mu2e::CaloClusterMC> CaloClusterMCCollection;

} 

#endif
