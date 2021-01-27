#include "MCDataProducts/inc/CaloClusterMC.hh"
#include "MCDataProducts/inc/CaloEDepMC.hh"
#include <numeric>

namespace mu2e {


   std::vector<CaloEDepMC> CaloClusterMC::energyDeposits() const
   {
       std::vector<CaloEDepMC> edeps;
       for (const auto& hit : hits_)
       {
          for (const auto& eDepInDigi : hit->energyDeposits())
          {
              auto it = edeps.begin();
              while (it != edeps.end()) {if (it->sim() == eDepInDigi.sim()) break;  ++it;}

              if (it!= edeps.end()) 
              {
                  it->addEDep(eDepInDigi.energyDep());
                  it->addEDepG4(eDepInDigi.energyDepG4());
                  it->addTime(eDepInDigi.time());
                  it->addMom(eDepInDigi.momentumIn());
              }
              else 
              {
                  edeps.emplace_back(CaloEDepMC(eDepInDigi.sim(),eDepInDigi.energyDep(),eDepInDigi.energyDepG4(),
                                                eDepInDigi.time(),eDepInDigi.momentumIn(),eDepInDigi.rel()));                
              }
          }
       } 
       std::sort(edeps.begin(),edeps.end(),[](const auto& a, const auto& b){return a.energyDep() > b.energyDep();});
       return edeps;
   }   


   float CaloClusterMC::totalEnergyDep() const 
   {
      auto sumEdep = [](float sum, const auto& hit){return sum +=hit->totalEnergyDep();};
      return std::accumulate(hits_.begin(),hits_.end(),0.0f,sumEdep);
   }

   float CaloClusterMC::totalEnergyDepG4() const 
   {
      auto sumEdepG4 = [](float sum, const auto& hit){return sum +=hit->totalEnergyDepG4();};
      return std::accumulate(hits_.begin(),hits_.end(),0.0f,sumEdepG4);
   }

}
