//
// Free function to construct version 5 of the TTracker
//
// Original author David N. Brown (Louisville), Nov 2017, based on
// previous versions by KLG and RKK.

//
// Notes
//
// 1)  The v5 in this function name says that this is the third way we
//     have implemented a single TTracker design in G4.  It does not refer
//     to alternate designs of the TTracker.
//
//     This version makes logical mother volumes per plane and per
//     panel and places panels in plane and straws in panel
//     This version places individual panels and planes, in order to 
//     work with alignment.
//
// 2) This function can build the TTracker designs described in:
//      Mu2eG4/test/ttracker_v2.txt   - Adjust spacings to match Mu2e-doc-888-v2.
//      and beyond
//

// C++ includes
#include <iomanip>
#include <iostream>
#include <string>

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "G4Helper/inc/G4Helper.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/constructTTracker.hh"
#include "Mu2eG4/inc/ConstructTTrackerTDR.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"

// G4 includes
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4IntersectionSolid.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"


using namespace std;

namespace mu2e{

  VolumeInfo constructTTrackerv5( VolumeInfo const& ds3Vac,
                                  SimpleConfig const& config ){

    // Master geometry for the TTracker.
    //    TTracker const & ttracker = *(GeomHandle<TTracker>());

    // The more detailed version has its own function.
    ConstructTTrackerTDR tt(ds3Vac, config);
    return tt.motherInfo();
    // Temporary until I can do more...

  } // end of constructTTrackerv5

} // end namespace mu2e
