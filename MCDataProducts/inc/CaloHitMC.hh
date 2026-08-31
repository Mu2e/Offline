//
// This calss contains CaloMCEdep for each CaloDigi, matching deposited energy at MC level to reco digi hits
// Notes:
// CaloEdepMC: the list of CaloEDepMC is ordered by energy - most energetic first.
// time:       taken as the time of the most energetic contribution.
//
//
#ifndef MCDataProducts_CaloHitMC_hh
#define MCDataProducts_CaloHitMC_hh

#include "Offline/MCDataProducts/inc/CaloEDepMC.hh"
#include "Offline/MCDataProducts/inc/SimParticleRemapping.hh"
#include <vector>

namespace mu2e
{

   class CaloHitMC
   {
       public:
          CaloHitMC() : edeps_(), crystalID_(-1) {};
          CaloHitMC(const std::vector<CaloEDepMC>& edeps, int crystalId) :
           edeps_(edeps),
           crystalID_(crystalId)
          {};
          CaloHitMC(const std::vector<CaloEDepMC>&& edeps, int crystalId) :
           edeps_(std::move(edeps)),
           crystalID_(crystalId)
          {};

          void resetSim(SimParticleRemapping const& remap);

          const std::vector<CaloEDepMC>& energyDeposits  ()           const {return edeps_;       }
          std::vector<CaloEDepMC>&       energyDeposits  ()                 {return edeps_;       }
          const CaloEDepMC&              energyDeposit   (unsigned i) const {return edeps_.at(i); }
          unsigned                       nParticles      ()           const {return edeps_.size();}
          float                          time            ()           const {return edeps_.empty() ? 0.0 : edeps_.at(0).time();}
          float                          totalEnergyDep  ()           const;
          float                          totalEnergyDepG4()           const;
          int                            crystalID       ()           const { return crystalID_;}

       private:
          std::vector<CaloEDepMC> edeps_;
          int                     crystalID_;
   };

   using  CaloHitMCCollection = std::vector<mu2e::CaloHitMC>;
}

#endif

