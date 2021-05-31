// Free function to create world mother volume and partly fill it with
// dirt around the formal hall box.
//
//
// Original author KLG based on Mu2eWorld constructDirt
// Updated by Andrei Gaponenko.

// art includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// Mu2e includes.
#include "Mu2eG4/inc/constructWorldVolume.hh"
#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4Box.hh"

using namespace std;

namespace mu2e {

  VolumeInfo constructWorldVolume(const SimpleConfig &config) {
    // A helper class.
    MaterialFinder materialFinder(config);

    // Dimensions and material of the world.
    G4Material* worldMaterial = materialFinder.get("world.materialName");
    G4Material* dirtMaterial  = materialFinder.get("dirt.overburdenMaterialName");

    G4GeometryOptions * geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "World"    , "world.box" );
    geomOptions->loadEntry( config, "worldDirt", "world.dirt");

    GeomHandle<WorldG4> world;

    VolumeInfo worldInfo = nestBox("World", world->halfLengths(),
                                   worldMaterial, 0, G4ThreeVector(),
                                   nullptr,
                                   0, G4Colour::Red() );

    //----------------------------------------------------------------
    // Here we create dirt slabs and place them in the World volume
    // around the formal "hall" box.

    // Dirt slab at the bottom
    {
      std::vector<double> hs(world->halfLengths());
      // world coordinates
      const double dirtYmin = -world->halfLengths()[1];
      const double dirtYmax = -world->hallFormalHalfSize()[1] + world->hallFormalCenterInWorld().y();
      hs[1] = (dirtYmax - dirtYmin)/2;

      CLHEP::Hep3Vector centerInWorld(3);
      centerInWorld[0] = 0.;
      centerInWorld[1] = (dirtYmax + dirtYmin)/2;
      centerInWorld[2] = 0;

      nestBox("worldDirtBottom", hs, dirtMaterial, 0, centerInWorld,
              worldInfo,
              0, G4Color::Magenta(), "worldDirt" );
    }

    // The height parameters are common to all 4 side slabs
    const double dirtYmin = -world->hallFormalHalfSize()[1] + world->hallFormalCenterInWorld().y();
    const double dirtYmax =  world->dirtG4Ymax();
    const double dirtCenterY = (dirtYmax + dirtYmin)/2;
    const double dirtHalfSizeY = (dirtYmax - dirtYmin)/2;

    // NW: slab covering the North corner and extending to the West
    {
      const double xmax = world->halfLengths()[0];
      const double xmin = world->hallFormalCenterInWorld().x() - world->hallFormalHalfSize()[0];
      const double zmin = -world->halfLengths()[2];
      const double zmax = world->hallFormalCenterInWorld().z() - world->hallFormalHalfSize()[2];
      std::vector<double> hs(3);
      hs[0] = (xmax - xmin)/2;
      hs[1] = dirtHalfSizeY;
      hs[2] = (zmax - zmin)/2;

      nestBox("worldDirtNW", hs, dirtMaterial, 0,
              CLHEP::Hep3Vector((xmax+xmin)/2, dirtCenterY, (zmax+zmin)/2),
              worldInfo,
              0, G4Color::Magenta(), "worldDirt" );
    }

    // SW: slab covering the West corner and extending to the South
    {
      const double xmin = -world->halfLengths()[0];
      const double xmax =  world->hallFormalCenterInWorld().x() - world->hallFormalHalfSize()[0];
      const double zmin = -world->halfLengths()[2];
      const double zmax = world->hallFormalCenterInWorld().z() + world->hallFormalHalfSize()[2];
      std::vector<double> hs(3);
      hs[0] = (xmax - xmin)/2;
      hs[1] = dirtHalfSizeY;
      hs[2] = (zmax - zmin)/2;

      nestBox("worldDirtSW", hs, dirtMaterial, 0,
              CLHEP::Hep3Vector((xmax+xmin)/2, dirtCenterY, (zmax+zmin)/2),
              worldInfo,
              0, G4Color::Magenta(), "worldDirt" );
    }

    // SE: slab covering the South corner and extending to the East
    {
      const double xmin = -world->halfLengths()[0];
      const double xmax =  world->hallFormalCenterInWorld().x() + world->hallFormalHalfSize()[0];
      const double zmax = +world->halfLengths()[2];
      const double zmin =  world->hallFormalCenterInWorld().z() + world->hallFormalHalfSize()[2];
      std::vector<double> hs(3);
      hs[0] = (xmax - xmin)/2;
      hs[1] = dirtHalfSizeY;
      hs[2] = (zmax - zmin)/2;

      nestBox("worldDirtSE", hs, dirtMaterial, 0,
              CLHEP::Hep3Vector((xmax+xmin)/2, dirtCenterY, (zmax+zmin)/2),
              worldInfo,
              0, G4Color::Magenta(), "worldDirt" );
    }

    // NE: slab covering the East corner and extending to the North
    {
      const double xmax = +world->halfLengths()[0];
      const double xmin =  world->hallFormalCenterInWorld().x() + world->hallFormalHalfSize()[0];
      const double zmax = +world->halfLengths()[2];
      const double zmin =  world->hallFormalCenterInWorld().z() - world->hallFormalHalfSize()[2];
      std::vector<double> hs(3);
      hs[0] = (xmax - xmin)/2;
      hs[1] = dirtHalfSizeY;
      hs[2] = (zmax - zmin)/2;

      nestBox("worldDirtNE", hs, dirtMaterial, 0,
              CLHEP::Hep3Vector((xmax+xmin)/2, dirtCenterY, (zmax+zmin)/2),
              worldInfo,
              0, G4Color::Magenta(), "worldDirt" );
    }

    //----------------------------------------------------------------
    return worldInfo;

  } // constructWorldVolume()

} // namespace mu2e
