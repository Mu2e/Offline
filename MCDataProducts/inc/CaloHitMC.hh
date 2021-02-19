//
// This calss contains CaloMCEdep for each CaloDigi, matching deposited energy at MC level to reco digi hits
// Notes:
// CaloEdepMC: the list of CaloEDepMC is ordered by energy - most energetic first. 
// time:       taken as the time of the most energetic contribution.
//
//
#ifndef MCDataProducts_CaloHitMC_hh
#define MCDataProducts_CaloHitMC_hh

#include "MCDataProducts/inc/CaloEDepMC.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"
#include <vector>

namespace mu2e
{   
   
   class CaloHitMC
   {
       public:                    
          CaloHitMC()                                     : edeps_()      {};
          CaloHitMC(const std::vector<CaloEDepMC>& edeps) : edeps_(edeps) {};
	  // reset SimParticle Ptrs
          void resetSim(SimParticleRemapping const& remap) {
	    for(auto& edep : edeps_){
	      auto newsim = remap.at(edep.sim());
	      edep.resetSim(newsim);
	    }
	  }
        
          const std::vector<CaloEDepMC>& energyDeposits  ()           const {return edeps_;       }
          std::vector<CaloEDepMC>& energyDeposits  ()           {return edeps_;       }
          const CaloEDepMC&              energyDeposit   (unsigned i) const {return edeps_.at(i); } 
          unsigned                       nParticles      ()           const {return edeps_.size();}
          float                          time            ()           const {return edeps_.at(0).time();}
          float                          totalEnergyDep  ()           const;
          float                          totalEnergyDepG4()           const;

       private:
          std::vector<CaloEDepMC> edeps_; 
   };  
   
   using  CaloHitMCCollection = std::vector<mu2e::CaloHitMC>;
} 

#endif

