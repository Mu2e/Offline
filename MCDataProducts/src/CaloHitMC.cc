#include "Offline/MCDataProducts/inc/CaloHitMC.hh"
#include <numeric>

namespace mu2e {

   float CaloHitMC::totalEnergyDep() const
   {
      auto sumEdep = [](float sum, const CaloEDepMC& edep){return sum +=edep.energyDep();};
      return std::accumulate(edeps_.begin(),edeps_.end(),0.0f,sumEdep);
   }

   float CaloHitMC::totalEnergyDepG4() const
   {
      auto sumEdepG4 = [](float sum, const CaloEDepMC& edep){return sum +=edep.energyDepG4();};
      return std::accumulate(edeps_.begin(),edeps_.end(),0.0f,sumEdepG4);
   }

   void CaloHitMC::resetSim(const SimParticleRemapping& remap)
   {
      for (auto& edep : edeps_)
      {
         auto newsim = remap.at(edep.sim());
         edep.resetSim(newsim);
      }
   }
}
