//
// Free function to create a geant4 test environment geometry
// N equally spaced boxes resembling in spacing the Mu2e stopping target
//
// Original author KLG
//
// Notes:
//
//

#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// Mu2e includes.

#include "Offline/Mu2eG4/inc/MaterialFinder.hh"
#include "Offline/Mu2eG4/inc/InitEnvToolBase.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"

// G4 includes
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "CLHEP/Units/SystemOfUnits.h"

// using namespace std;

namespace mu2e {

  class ConstructEnvNBoxes: public InitEnvToolBase {
  public:
    ConstructEnvNBoxes(const fhicl::ParameterSet& PSet);
    int construct(VolumeInfo const& ParentVInfo, SimpleConfig const& Config);
  };


//-----------------------------------------------------------------------------
  ConstructEnvNBoxes::ConstructEnvNBoxes(const fhicl::ParameterSet& PSet) {
    _name = "NBoxes";
  }
//-----------------------------------------------------------------------------
  int ConstructEnvNBoxes::construct(VolumeInfo const& parentVInfo, SimpleConfig const& _config) {

    const G4bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible");
    const G4bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck");
    const G4bool placePV             = true;

    // Extract tube information from the config file.

    G4bool boxVisible        = _config.getBool("box.visible",true);
    G4bool boxSolid          = _config.getBool("box.solid",true);

    std::vector<double> boxParams;
    _config.getVectorDouble( "box.halfLengths", boxParams);

    MaterialFinder materialFinder(_config);
    G4Material* boxMaterial = materialFinder.get("box.materialName");

    const G4ThreeVector firstBoxCenterInWorld(_config.getHep3Vector("box.centerInWorld"));

    G4Colour  orange  (.75, .55, .0);

    std::ostringstream vsBoxNumber("");
    std::string vBoxName;

    G4int startingBoxCopyNumber = _config.getInt("box.copyNumber",2);

    G4int        numberOfBoxes  = _config.getInt("box.numberOfBoxes");
    G4double     boxSpacing     = _config.getDouble("box.spacing");

    std::string mNamePrefix("Box");

    for( G4int nb = 0; nb<numberOfBoxes; ++nb ) {

      G4int boxCopyNumber = startingBoxCopyNumber+nb;
      vsBoxNumber.str("");
      vsBoxNumber.width(2);
      vsBoxNumber.fill('0');
      vsBoxNumber << boxCopyNumber;
      vBoxName = mNamePrefix + vsBoxNumber.str();

      G4ThreeVector boxCenterInWorld = firstBoxCenterInWorld +
        G4ThreeVector(0.,0.,nb*boxSpacing); // we space them along the z axis

      VolumeInfo boxVInfo(nestBox( vBoxName,
                                   boxParams,
                                   boxMaterial,
                                   nullptr, // no rotation
                                   boxCenterInWorld,
                                   parentVInfo,
                                   boxCopyNumber,
                                   // we assign a non 0 copy nuber for
                                   // volume tracking purposes
                                   boxVisible,
                                   orange,
                                   boxSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ));

    }

    return 0;
  }

}

DEFINE_ART_CLASS_TOOL(mu2e::ConstructEnvNBoxes)
