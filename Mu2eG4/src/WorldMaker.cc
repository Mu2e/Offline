//
// The Mu2e version of G4VUserDetectorConstruction.
//
//
// Original author Rob Kutschke
//
// This receives calls from G4 and forwards them to
// the class Mu2eUniverse that actually constructs the
// Mu2e world.  This class manages the lifetime of
// the Mu2eUniverse and ConstructMaterials objects.
//
// Modified by Lisa Goodenough 7/11/17 to add MT functionality
// through G4 ConstructSDandField method, which instantiates the
// thread-local instances of SensitiveDetectors in the threads
//

#include "Mu2eG4/inc/WorldMaker.hh"

#include "Mu2eG4/inc/Mu2eStudyWorld.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include <iostream>


namespace mu2e {

  // Construct() is called by G4 and just returns world physical volume
  template <typename WorldType, typename MaterialsType>
  G4VPhysicalVolume* WorldMaker<WorldType, MaterialsType>::Construct(){
    // Clean old geometry, if any
    Clean();

    _materials->construct();
    return _world->construct();
  }

  // ConstructSDandField() is called by G4 and instantiates the SensitiveDetectors
  template <typename WorldType, typename MaterialsType>
  void WorldMaker<WorldType, MaterialsType>::ConstructSDandField(){

    _world->constructSDandField();
  }


  // Clean old geometry, if any.
  template <typename WorldType, typename MaterialsType>
  void WorldMaker<WorldType, MaterialsType>::Clean(){

    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
  }


  template class WorldMaker<Mu2eStudyWorld,ConstructMaterials>;
  template class WorldMaker<Mu2eWorld,ConstructMaterials>;

} // end namespace mu2e
