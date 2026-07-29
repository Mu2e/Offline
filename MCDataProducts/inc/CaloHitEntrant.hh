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
// Grouping deposits by entrant recovers true shower membership for
// clustering truth definitions. Analysis-level choices (purity cuts,
// ambiguity handling, cluster ID assignment) are deliberately left to
// consumers; this product records only the ancestry facts.
//
// The collection is index-parallel to the CaloHitMCCollection it was
// produced from.
//

#include "canvas/Persistency/Common/Ptr.h"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include <vector>

namespace mu2e
{
   struct CaloHitEntrant
   {
      std::vector<art::Ptr<SimParticle>> entrants;
   };

   using CaloHitEntrantCollection = std::vector<mu2e::CaloHitEntrant>;
}

#endif
