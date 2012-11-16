#ifndef Mu2eG4_WorldMaker_hh
#define Mu2eG4_WorldMaker_hh
//
// The Mu2e version of G4VUserDetectorConstruction.
//
// $Id: WorldMaker.hh,v 1.5 2012/11/16 23:45:06 genser Exp $
// $Author: genser $
// $Date: 2012/11/16 23:45:06 $
//
// Original author Rob Kutschke
//
// This receives calls from G4 and forwards them to
// the class Mu2eUniverse that actually constructs the
// Mu2e world.  This class manages the lifetime of
// the Mu2eUniverse and ConstructMaterials objects.
//

//#include <string>
#include <memory>

// Mu2e includes
#include "Mu2eG4/inc/ConstructMaterials.hh"

// Forward references.
// class G4VPhysicalVolume;

// G4 includes
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VUserDetectorConstruction.hh"

namespace mu2e {

  template <typename T> class WorldMaker : public G4VUserDetectorConstruction
  {
  public:

    WorldMaker():
      _materials(),
      _world(){
    }
    ~WorldMaker(){}

    // This is the required method prescribed by G4.
    G4VPhysicalVolume* Construct(){

      // Clean old geometry, if any
      Clean();

      _materials = std::auto_ptr<ConstructMaterials>(new ConstructMaterials());
      _world     = std::auto_ptr<T>(new T());

      _materials->construct();

      return _world ->construct();
    }

    // Accessors.
    T const* getWorld()     { return _world.get(); }

  private:


    // Clean old geometry, if any.
    void Clean(){

      _materials.release();
      _world.release();

      G4GeometryManager::GetInstance()->OpenGeometry();
      G4PhysicalVolumeStore::GetInstance()->Clean();
      G4LogicalVolumeStore::GetInstance()->Clean();
      G4SolidStore::GetInstance()->Clean();

    }

    std::auto_ptr<ConstructMaterials> _materials;
    std::auto_ptr<T>                  _world;

  };

} // end namespace mu2e
#endif /* Mu2eG4_WorldMaker_hh */

