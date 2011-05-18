//
// Free function to create world mother volume
//
// $Id: constructWorldVolume.cc,v 1.4 2011/05/18 14:21:44 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/18 14:21:44 $
//
// Original author KLG based on Mu2eWorld constructDirt
//
// Notes:
// Construct the (top) World Mother Volume

// Mu2e includes.
#include "Mu2eG4/inc/constructWorldVolume.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Box.hh"

using namespace std;

namespace mu2e {

  VolumeInfo constructWorldVolume( SimpleConfig const * const _config
                                   ){

    // A helper class.
    MaterialFinder materialFinder(*_config);

    // Dimensions and material of the world.
    vector<double> worldHLen;
    _config->getVectorDouble("world.halfLengths", worldHLen, 3);
    G4Material* worldMaterial = materialFinder.get("world.materialName");

    bool worldBoxVisible = _config->getBool("world.boxVisible",true);
    bool worldBoxSolid   = _config->getBool("world.boxSolid",false);

    bool const forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool const placePV             = true;

    VolumeInfo worldInfo;

    worldInfo.name    = "World";

    worldInfo.solid   = new G4Box( worldInfo.name, worldHLen[0], worldHLen[1], worldHLen[2] );

    finishNesting(worldInfo,
                  worldMaterial,
                  0,
                  G4ThreeVector(),
                  0,
                  0,
                  worldBoxVisible,
                  G4Colour::Red(),
                  worldBoxSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  false
                  );

    return worldInfo;

  }

}
