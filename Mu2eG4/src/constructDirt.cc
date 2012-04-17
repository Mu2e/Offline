// Constructs dirt overburden inside the formal hall box.
// Note that there is also dirt around the hall box.
//
// $Id: constructDirt.cc,v 1.9 2012/04/17 19:56:56 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/17 19:56:56 $
//
// Original author KLG based on Mu2eWorld constructDirt
// Updated by Andrei Gaponenko.

// Mu2e includes.
#include "Mu2eG4/inc/constructDirt.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Paraboloid.hh"

using namespace std;

namespace mu2e {

  void constructDirt(const VolumeInfo& parent, const SimpleConfig& config) {

    // Here we can e.g. place dirt on top of the hall ceiling.


  } // constructDirt()

} // namespace mu2e
