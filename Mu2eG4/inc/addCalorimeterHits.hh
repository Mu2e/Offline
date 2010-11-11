#ifndef AddCalorimeterHits_HH
#define AddCalorimeterHits_HH
//
// Populate output collection for calorimeter
//
// $Id: addCalorimeterHits.hh,v 1.2 2010/11/11 21:19:19 genser Exp $
// $Author: genser $
// $Date: 2010/11/11 21:19:19 $
//
// Original author Ivan Logashenko
//

#include "ToyDP/inc/CaloHitCollection.hh"
#include "ToyDP/inc/CaloHitMCTruthCollection.hh"
#include "ToyDP/inc/CaloCrystalHitMCTruthCollection.hh"

class G4Event;

namespace mu2e{

  // Public entry point.
  void addCalorimeterHits ( const G4Event * event, 
			    CaloHitCollection& caloHits, 
			    CaloHitMCTruthCollection& caloHitsMCTruth,
                            CaloCrystalHitMCTruthCollection& caloCrystalHitsMCTruth);

}  // end namespace mu2e

#endif


