#include "MCDataProducts/inc/CaloClusterMC.hh"
#include <numeric>

namespace mu2e {

   bool CaloClusterMC::isConversion() const 
   {
      for (auto& edep : edeps_)
      {
         auto parent(edep.sim());
         while (parent->hasParent()) parent = parent->parent();                     
	 if (parent->genParticle() && parent->genParticle()->generatorId().isConversion() ) return true;
      }    		
      return false;
   }


   bool CaloClusterMC::isOnlyConversion() const 
   {
      for (auto& edep : edeps_)
      {
         auto parent(edep.sim());
         while (parent->hasParent()) parent = parent->parent();                     
         if ( !(parent->genParticle() && parent->genParticle()->generatorId().isConversion()) ) return false;
      }    		
      return false;
   }


   float CaloClusterMC::totalEnergyDep() const 
   {
      auto sumEdep = [](float sum, const CaloEDepMC& edep){return sum +=edep.energyDep();};
      return std::accumulate(edeps_.begin(),edeps_.end(),0.0f,sumEdep);
   }


   float CaloClusterMC::totalEnergyDepG4() const 
   {
      auto sumEdepG4 = [](float sum, const CaloEDepMC& edep){return sum +=edep.energyDepG4();};
      return std::accumulate(edeps_.begin(),edeps_.end(),0.0f,sumEdepG4);
   }

}
