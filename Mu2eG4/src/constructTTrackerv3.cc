//
// Free function to construct version 3 of the TTracker
//
// $Id: constructTTrackerv3.cc,v 1.38 2014/04/11 04:42:06 genser Exp $
// $Author: genser $
// $Date: 2014/04/11 04:42:06 $
//
// Original author KLG based on RKK's version using different methodology
//
// Notes
//
// 1)  The v3 in this function name says that this is the third way we
//     have implemented a single TTracker design in G4.  It does not refer
//     to alternate designs of the TTracker.
//
//     This version makes logical mother volumes per plane and per
//     panel and places panels in plane and straws in panel
//     It has only one panel/plane logical volume placed several times
//     This version has a negligeable construction time and a much smaler memory footprint
//
// 2) This function can build the TTracker designs described in:
//      Mu2eG4/test/ttracker_meco.txt - The MECO design, uniform plane spacing
//      Mu2eG4/test/ttracker_v0.txt   - The first Aseet version, pairs of planes form stations
//                                      but one layer of straws per panel (called a panel in this code)
//      Mu2eG4/test/ttracker_v1.txt   - v0 but with with two layers of straws per panel
//      Mu2eG4/test/ttracker_v2.txt   - Adjust spacings to match Mu2e-doc-888-v2.
//
// 3) This function does not know how to build the TTracker described in:
//       Mu2eG4/test/ttracker_v3.txt - Detail support model and detailed layering of straws
//    This geometry can be detected by the method by
//
//      if ( ttracker.getSupportModel() == SupportModel::detailedv0 ) ....
//
//    If this geometry is detected, this function call through to constructTTrackerv3Detailed.cc

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

  VolumeInfo constructTTrackerv3( VolumeInfo const& ds3Vac,
                                  SimpleConfig const& config ){

    // Master geometry for the TTracker.
    TTracker const & ttracker = *(GeomHandle<TTracker>());

    // The more detailed version has its own function.
    if ( ttracker.getSupportModel() == SupportModel::detailedv0 ) {
      ConstructTTrackerTDR tt(ds3Vac, config);
      return tt.motherInfo();
    }
      
    throw cet::exception("GEOM")  << "Unsupported Tracker Config\n";

  } // end of constructTTrackerv3

} // end namespace mu2e
