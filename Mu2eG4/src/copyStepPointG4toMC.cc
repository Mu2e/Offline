//
// Helper function to copy StepPointG4s to output StepPointMC collection
//
// Original author Ivan Logashenko
//

// Mu2e includes
#include "Mu2eG4/inc/StepPointG4.hh"
#include "Mu2eG4/inc/copyStepPointG4toMC.hh"

// G4 includes
#include "G4Event.hh"
#include "G4SDManager.hh"

namespace mu2e{

  void  copyStepPointG4toMC ( const G4Event* g4event,
                              const std::string name,
                              StepPointMCCollection& outputHits  ){

    // G4 Hit collections for this event.

    G4HCofThisEvent* hce = g4event->GetHCofThisEvent();

    // Get the collection ID for the VD hits.

    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4int colId = SDman->GetCollectionID(name);

    // Copy hits from G4 collection to output Mu2e collection

    if ( colId >= 0 && hce != 0 ){

      StepPointG4Collection* hits = static_cast<StepPointG4Collection*>(hce->GetHC(colId));
      G4int nHits = hits->entries();

      for (G4int i=0;i<nHits;i++) {
        StepPointG4* h = (*hits)[i];
        outputHits.push_back( h->hit() );
      }

    }

  }

} // end namespace mu2e
