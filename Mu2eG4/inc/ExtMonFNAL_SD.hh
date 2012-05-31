#ifndef Mu2eG4_ExtMonFNAL_SD_hh
#define Mu2eG4_ExtMonFNAL_SD_hh
//
// $Id: ExtMonFNAL_SD.hh,v 1.4 2012/05/31 17:08:02 genser Exp $
// $Author: genser $
// $Date: 2012/05/31 17:08:02 $
//
// Define a sensitive detector for the FNAL extinction monitor
// Derived from the CaloCrystalSD code.
// Andrei Gaponenko, 2011
//

#include <map>
#include <vector>

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"

// Art includes
#include "art/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

namespace mu2e {

  class ExtMonFNAL_SD : public Mu2eSensitiveDetector{

  public:

    ExtMonFNAL_SD(G4String const name, SimpleConfig const & config);

    G4bool ProcessHits(G4Step*, G4TouchableHistory*);

  };

} // namespace mu2e

#endif /* Mu2eG4_ExtMonFNAL_SD_hh */
