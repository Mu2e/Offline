//
// Free function to create a geant4 test environment geometry
//
// $Id: constructStudyEnv_v004.cc,v 1.1 2013/02/27 03:49:02 genser Exp $
// $Author: genser $
// $Date: 2013/02/27 03:49:02 $
//
// Original author KLG 
//
// Notes:
//
// one can nest volume inside other volumes if needed
// see other construct... functions for examples
//

// Mu2e includes.

#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/constructStudyEnv_v004.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

using namespace std;

namespace mu2e {

  void constructStudyEnv_v004( VolumeInfo   const & parentVInfo,
                               SimpleConfig const & _config
                               ){

    const bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible");
    const bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck");
    const bool placePV             = true;

    // Extract box information from the config file.

    G4bool boxVisible        = _config.getBool("box.visible",true);
    G4bool boxSolid          = _config.getBool("box.solid",true);

    vector<double> boxParams;
    _config.getVectorDouble( "box.halfLengths", boxParams);

    MaterialFinder materialFinder(_config);
    G4Material* boxMaterial = materialFinder.get("box.wallMaterialName");

    const G4ThreeVector boxCenterInWorld(_config.getHep3Vector("box.centerInWorld"));

    G4Colour  orange  (.75, .55, .0);

    VolumeInfo boxVInfo(nestBox( "Box",
                                 boxParams,
                                 boxMaterial,
                                 0, // no rotation
                                 boxCenterInWorld,
                                 parentVInfo,
                                 _config.getInt("box.copyNumber",2), 
                                 // we assign a non 0 copy nuber for
                                 // volume tracking purposes
                                 boxVisible,
                                 orange,
                                 boxSolid,
                                 forceAuxEdgeVisible,
                                 placePV,
                                 doSurfaceCheck
                                 ));

    // Extract tube information from the config file.

    G4bool tubeVisible        = _config.getBool("tube.visible",true);
    G4bool tubeSolid          = _config.getBool("tube.solid",true);


    TubsParams tubeParams( _config.getDouble("tube.rIn"),
                           _config.getDouble("tube.rOut"),
                           _config.getDouble("tube.halfLength"));

    G4Material* tubeMaterial = materialFinder.get("tube.wallMaterialName");

    const G4ThreeVector tubeCenterInWorld(_config.getHep3Vector("tube.centerInWorld"));

    VolumeInfo tubeVInfo(nestTubs( "Tube",
                                   tubeParams,
                                   tubeMaterial,
                                   0, // no rotation
                                   tubeCenterInWorld,
                                   parentVInfo,
                                   _config.getInt("tube.copyNumber",3), 
                                   // we assign a non 0 copy nuber for
                                   // volume tracking purposes
                                   tubeVisible,
                                   orange,
                                   tubeSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ));






  } // constructStudyEnv_v004;

}
