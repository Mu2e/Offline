#include "MCDataProducts/inc/CaloClusterNewMC.hh"
#include <numeric>

namespace mu2e {

   bool CaloClusterNewMC::isConversion() const 
   {
      for (auto& edep : edeps_)
      {
         auto parent(edep.sim());
         while (parent->hasParent()) parent = parent->parent();                     
	 if (parent->genParticle() && parent->genParticle()->generatorId().isConversion() ) return true;
      }    		
      return false;
   }


   bool CaloClusterNewMC::isOnlyConversion() const 
   {
      for (auto& edep : edeps_)
      {
         auto parent(edep.sim());
         while (parent->hasParent()) parent = parent->parent();                     
         if ( !(parent->genParticle() && parent->genParticle()->generatorId().isConversion()) ) return false;
      }    		
      return false;
   }


   float CaloClusterNewMC::totalEnergyDep() const 
   {
      auto sumEdep = [](float sum, const CaloEDepMC& edep){return sum +=edep.eDep();};
      return std::accumulate(edeps_.begin(),edeps_.end(),0.0f,sumEdep);
   }


   float CaloClusterNewMC::totalEnergyDepG4() const 
   {
      auto sumEdepG4 = [](float sum, const CaloEDepMC& edep){return sum +=edep.eDepG4();};
      return std::accumulate(edeps_.begin(),edeps_.end(),0.0f,sumEdepG4);
   }

}
