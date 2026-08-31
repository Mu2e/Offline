#ifndef MCDataProducts_CaloHitEntrant_hh
#define MCDataProducts_CaloHitEntrant_hh
//
// Calo-entrant truth assignment for one CaloHitMC.
//
// For each MC energy deposit in the hit (aligned with
// CaloHitMC::energyDeposits()), stores the SimParticle that originated
// the shower on this disk: the highest ancestor in the Geant4 parent
// chain that also deposited energy in the same calorimeter disk.
// A particle with no such ancestor (e.g. a cross-disk secondary or a
// pileup particle arriving from outside) is its own calo-entrant.
//
// caloHitMC() points back to the hit this assignment was computed for,
// so the product is self-describing: consumers join through the Ptr
// (or verify it) instead of trusting that their configured
// CaloHitMCCollection is the one the producer read. The collection is
// also index-parallel to that CaloHitMCCollection by construction.
//
// Grouping deposits by entrant recovers true shower membership for
// clustering truth definitions. Analysis-level choices (purity cuts,
// ambiguity handling, cluster ID assignment) are deliberately left to
// consumers; this product records only the ancestry facts.
//

#include "canvas/Persistency/Common/Ptr.h"
#include "Offline/MCDataProducts/inc/CaloHitMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include <utility>
#include <vector>

namespace mu2e
{
   class CaloHitEntrant
   {
      public:
         CaloHitEntrant() = default;
         CaloHitEntrant(art::Ptr<CaloHitMC> caloHitMC,
                        std::vector<art::Ptr<SimParticle>> entrants) :
            caloHitMC_(std::move(caloHitMC)),
            entrants_ (std::move(entrants))
         {}

         const art::Ptr<CaloHitMC>&                caloHitMC() const {return caloHitMC_;}
         const std::vector<art::Ptr<SimParticle>>& entrants () const {return entrants_;}

      private:
         art::Ptr<CaloHitMC>                caloHitMC_;
         std::vector<art::Ptr<SimParticle>> entrants_;
   };

   using CaloHitEntrantCollection = std::vector<mu2e::CaloHitEntrant>;
}

#endif
