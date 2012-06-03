#ifndef Mu2eG4_ExtMonUCI_TofSD_hh
#define Mu2eG4_ExtMonUCI_TofSD_hh
//
// Define an extinction monitor TOF detector
//
// $Id: ExtMonUCITofSD.hh,v 1.2 2012/06/03 06:54:57 youzy Exp $
// $Author: youzy $
// $Date: 2012/06/03 06:54:57 $
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

  class ExtMonUCITofSD : public Mu2eSensitiveDetector{

  public:

    ExtMonUCITofSD(G4String const name, SimpleConfig const & config);

    G4bool ProcessHits(G4Step*, G4TouchableHistory*);

  };

} // namespace mu2e

#endif /* Mu2eG4_ExtMonUCI_TofSD_hh */
