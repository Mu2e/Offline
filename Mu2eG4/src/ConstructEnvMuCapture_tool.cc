//
// Free function to create a geant4 test environment geometry
//
// $Id: constructStudyEnv_v001.cc,v 1.3 2013/02/27 03:49:40 genser Exp $
// $Author: genser $
// $Date: 2013/02/27 03:49:40 $
//
// Original author KLG 
//
// Notes:
//
// one can nest volume inside other volumes if needed
// see other construct... functions for examples
//

#include "art/Utilities/ToolMacros.h"

// Mu2e includes.

#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/InitEnvToolBase.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

namespace mu2e {

  class ConstructEnvMuCapture: public InitEnvToolBase {
  public:
    ConstructEnvMuCapture(const fhicl::ParameterSet& PSet);
    ~ConstructEnvMuCapture();

    int construct(VolumeInfo const& ParentVInfo, SimpleConfig const& Config);
  };


//-----------------------------------------------------------------------------
  ConstructEnvMuCapture::ConstructEnvMuCapture(const fhicl::ParameterSet& PSet) {
    _name = "MuCapture";
  }

//-----------------------------------------------------------------------------
  ConstructEnvMuCapture::~ConstructEnvMuCapture() {
    _name = "MuCapture";
  }

//-----------------------------------------------------------------------------
  int ConstructEnvMuCapture::construct(VolumeInfo const& parentVInfo, SimpleConfig const& _config) {

    const bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible");
    const bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck");
    const bool placePV             = true;

    // Extract tube information from the config file.

    G4bool tubeVisible        = _config.getBool("tube.visible",true);
    G4bool tubeSolid          = _config.getBool("tube.solid",true);


    TubsParams tubeParams( _config.getDouble("tube.rIn"),
                           _config.getDouble("tube.rOut"),
                           _config.getDouble("tube.halfLength"),
                           _config.getDouble("tube.phi0")*CLHEP::degree,
                           _config.getDouble("tube.phiSpan")*CLHEP::degree );

    MaterialFinder materialFinder(_config);
    G4Material* tubeMaterial = materialFinder.get("tube.wallMaterialName");

    const G4ThreeVector tubeCenterInWorld(_config.getHep3Vector("tube.centerInWorld"));

    G4Colour  orange  (.75, .55, .0);

    double sgn(_config.getDouble("tube.sign"));

    G4Helper& helper(*art::ServiceHandle<G4Helper>());
    AntiLeakRegistry& reg(helper.antiLeakRegistry());

    CLHEP::HepRotationZ rotZ(_config.getDouble("tube.phiRotZ")*CLHEP::degree);
    CLHEP::HepRotationY rotY(M_PI);
    G4RotationMatrix* rotation  = (sgn < 0 )?
      reg.add(G4RotationMatrix(rotZ)) :
      reg.add(G4RotationMatrix(rotZ*rotY));


    VolumeInfo tubeVInfo(nestTubs( "Tube",
                                   tubeParams,
                                   tubeMaterial,
                                   rotation, // rotation
                                   tubeCenterInWorld,
                                   parentVInfo,
                                   _config.getInt("tube.copyNumber",2), 
                                   // we assign a non 0 copy nuber for
                                   // volume tracking purposes
                                   tubeVisible,
                                   orange,
                                   tubeSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ));
    return 0;
  }

}

DEFINE_ART_CLASS_TOOL(mu2e::ConstructEnvMuCapture)
