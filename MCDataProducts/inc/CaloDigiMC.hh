//
// This calss contains CaloMCEdep for each CaloDigi, matching deposited energy at MC level to reco digi hits
// Notes:
// CaloEdepMC: the list of CaloEDepMC is ordered by energy - most energetic first. 
// time:       taken as the time of the most energetic contribution.
//
//
#ifndef MCDataProducts_CaloDigiMC_hh
#define MCDataProducts_CaloDigiMC_hh

#include "MCDataProducts/inc/CaloEDepMC.hh"
#include <vector>

namespace mu2e
{   
   
   class CaloDigiMC
   {
       public:                    
          CaloDigiMC()                                     : edeps_()      {};
          CaloDigiMC(const std::vector<CaloEDepMC>& edeps) : edeps_(edeps) {};
                
          const std::vector<CaloEDepMC>& energyDeposits  ()           const {return edeps_;       }
          const CaloEDepMC&              energyDeposit   (unsigned i) const {return edeps_.at(i); } 
          unsigned                       nParticles      ()           const {return edeps_.size();}
          float                          time            ()           const {return edeps_.at(0).time();}
          float                          totalEnergyDep  ()           const;
          float                          totalEnergyDepG4()           const;

       private:
          std::vector<CaloEDepMC> edeps_; 
   };  
   
   using  CaloDigiMCCollection = std::vector<mu2e::CaloDigiMC>;
} 

#endif

