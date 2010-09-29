#ifndef AddCalorimeterHits_HH
#define AddCalorimeterHits_HH
//
// Populate output collection for calorimeter
//
// Original author Ivan Logashenko
//

#include "ToyDP/inc/CaloHitCollection.hh"
#include "ToyDP/inc/CaloHitMCTruthCollection.hh"

class G4Event;

namespace mu2e{

  // Public entry point.
  void addCalorimeterHits ( const G4Event *, 
			    CaloHitCollection &hits, 
			    CaloHitMCTruthCollection& mchits );

}  // end namespace mu2e

#endif


