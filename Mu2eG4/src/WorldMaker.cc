//
// The Mu2e version of G4VUserDetectorConstruction.
// 
// $Id: WorldMaker.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//
// This receives calls from G4 and forwards them to
// the class Mu2eWorld that actually constructs the
// Mu2e world.  This class manages the lifetime of
// the Mu2eWorld and ConstructMaterials objects.
//

// C++ includes
#include <iostream>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eG4/inc/WorldMaker.hh"
#include "Mu2eG4/inc/ConstructMaterials.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

// G4 includes
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4PVPlacement.hh"
#include "globals.hh"

using namespace std;

namespace mu2e {

  WorldMaker::WorldMaker():
    _materials(),
    _world(){
  }
  
  WorldMaker::~WorldMaker(){
  }

  // This is the callback used by G4.
  G4VPhysicalVolume* WorldMaker::Construct(){

    // Clean old geometry, if any
    Clean();

    // Get access to the master geometry system and its run time config.
    edm::Service<GeometryService> geom;
    SimpleConfig const& config = geom->config();

    _materials = auto_ptr<ConstructMaterials>(new ConstructMaterials());
    _world     = auto_ptr<Mu2eWorld>(new Mu2eWorld());

    _materials->construct();

    _worldInfo = _world ->construct();

    return _worldInfo->worldPhys;
  }


  // Clean old geometry, if any.
  void WorldMaker::Clean(){

    _materials.release();
    _world.release();

    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

  }  
  
  
} // end namespace mu2e
