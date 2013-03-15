#ifndef Mu2eG4_WorldMaker_hh
#define Mu2eG4_WorldMaker_hh
//
// The Mu2e version of G4VUserDetectorConstruction.
//
// $Id: WorldMaker.hh,v 1.8 2013/03/15 20:11:23 russo Exp $
// $Author: russo $
// $Date: 2013/03/15 20:11:23 $
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

#include "boost/noncopyable.hpp"

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

  template <typename WorldType, typename MaterialsType=ConstructMaterials>
  class WorldMaker : public G4VUserDetectorConstruction,
                     private boost::noncopyable
  {
  public:

    explicit WorldMaker(std::unique_ptr<WorldType> pw = std::unique_ptr<WorldType>(new WorldType()),
                        std::unique_ptr<MaterialsType> pm = std::unique_ptr<MaterialsType>(new ConstructMaterials())) :
      _materials(std::move(pm)),
      _world(std::move(pw))
    {
    }
    ~WorldMaker(){}

    // This is the required method prescribed by G4.
    G4VPhysicalVolume* Construct(){

      // Clean old geometry, if any
      Clean();

      _materials->construct();

      return _world ->construct();
    }

    // Accessors.
    WorldType const* getWorld()     { return _world.get(); }

  private:

    // Clean old geometry, if any.
    void Clean(){

      G4GeometryManager::GetInstance()->OpenGeometry();
      G4PhysicalVolumeStore::GetInstance()->Clean();
      G4LogicalVolumeStore::GetInstance()->Clean();
      G4SolidStore::GetInstance()->Clean();

    }

    std::unique_ptr<MaterialsType> _materials;
    std::unique_ptr<WorldType>     _world;

  };

} // end namespace mu2e
#endif /* Mu2eG4_WorldMaker_hh */
