#ifndef Mu2eG4_WorldMaker_hh
#define Mu2eG4_WorldMaker_hh
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


#include <memory>

#include "boost/noncopyable.hpp"

// Mu2e includes
#include "Mu2eG4/inc/ConstructMaterials.hh"

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

    explicit WorldMaker(std::unique_ptr<WorldType> pw, std::unique_ptr<MaterialsType> pm)
      :
      _materials(std::move(pm)),
      _world(std::move(pw))
    {
    }

    ~WorldMaker(){}


    // These are required methods prescribed by G4.

    // Construct() is called by GEANT and just returns world physical volume
    // Construct() method should contain definition of materials, volumes and visualization attributes
    G4VPhysicalVolume* Construct();

    // Given sensitive detector class objects should be thread-local, instantiation of such
    // thread-localclasses should be implemented in this method ConstructSDandField()
    void ConstructSDandField();

    // Accessors
    WorldType const* getWorld()     { return _world.get(); }

  private:

    // Clean old geometry, if any.
    void Clean();
    std::unique_ptr<MaterialsType> _materials;
    std::unique_ptr<WorldType>     _world;

  };

} // end namespace mu2e
#endif /* Mu2eG4_WorldMaker_hh */
