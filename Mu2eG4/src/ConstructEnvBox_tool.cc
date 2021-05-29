//
// Free function to create a geant4 test environment geometry
//
//
// Original author KLG 
//
// Notes:
//
// one can nest volume inside other volumes if needed
// see other construct... functions for examples
//

#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// Mu2e includes.

#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/InitEnvToolBase.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// G4 includes
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

namespace mu2e {

  class ConstructEnvBox: public InitEnvToolBase {
  public:
    ConstructEnvBox(const fhicl::ParameterSet& PSet);
    ~ConstructEnvBox();

    int construct(VolumeInfo const& ParentVInfo, SimpleConfig const& Config);
  };


//-----------------------------------------------------------------------------
  ConstructEnvBox::ConstructEnvBox(const fhicl::ParameterSet& PSet) {
    _name = "Box";
  }

//-----------------------------------------------------------------------------
  ConstructEnvBox::~ConstructEnvBox() {
    _name = "Box";
  }

//-----------------------------------------------------------------------------
  int ConstructEnvBox::construct(VolumeInfo const& parentVInfo, SimpleConfig const& _config) {

    const bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible");
    const bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck");
    const bool placePV             = true;

    // Extract tube information from the config file.

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


    return 0;
  }

}

DEFINE_ART_CLASS_TOOL(mu2e::ConstructEnvBox)
