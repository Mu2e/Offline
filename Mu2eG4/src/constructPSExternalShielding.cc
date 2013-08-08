// David Norvil Brown, August 2013

#include "Mu2eG4/inc/constructPSExternalShielding.hh"

#include "CLHEP/Vector/TwoVector.h"

#include "ProductionSolenoidGeom/inc/PSExternalShielding.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
//#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Orb.hh"
#include "G4TwoVector.hh"
#include "CLHEP/Vector/Rotation.h"
#include "G4NistManager.hh"

#include <vector>
#include <algorithm>


namespace mu2e {

  //================================================================

  void constructPSExternalShielding(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<PSExternalShielding> psexs;

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    const bool placePV             = true;

    //----------------------------------------------------------------


//     VolumeInfo shield = nestBox("PSExternalShielding",
// 				      psexs->getBox(),
// 				      findMaterialOrThrow(psexs->materialName()),
// 				      0,
// 				      psexs->centerOfShield() - parent.centerInMu2e(),
// 				      parent,
// 				      0,
// 				      config.getBool("PSExternalShielding.visible"),
// 				      G4Colour::Magenta(),
// 				      config.getBool("PSExternalShielding.solid"),
// 				      forceAuxEdgeVisible,
// 				      placePV,
// 				      doSurfaceCheck
// 				      );

    // convert internally held std::vector of Hep2Vectors and convert to
    // G4TwoVectors for G4

    std::vector<G4TwoVector> verticalOutline; 
    std::copy(psexs->externalShieldOutline().begin(),
	      psexs->externalShieldOutline().end(),
	      std::back_inserter(verticalOutline));

    VolumeInfo shield("PSExternalShielding",
		      psexs->centerOfShield()-parent.centerInMu2e(),
		      parent.centerInWorld);

    shield.solid = new G4ExtrudedSolid(shield.name,
				       verticalOutline,
				       psexs->getLength()/2.,
				       G4TwoVector(0,0), 1.,
				       G4TwoVector(0,0), 1.);

    // no rotation, but need one for next call
    static CLHEP::HepRotation rotat(CLHEP::HepRotation::IDENTITY);

     finishNesting(shield,
 		  findMaterialOrThrow(psexs->materialName()),
 		  &rotat,
 		  shield.centerInParent,
 		  parent.logical,
 		  0,
 		  config.getBool("PSExternalShielding.visible"),
 		  G4Colour::Magenta(),
 		  config.getBool("PSExternalShielding.solid"),
 		  forceAuxEdgeVisible,
 		  placePV,
 		  doSurfaceCheck);

  }


}
