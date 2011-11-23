#include "Mu2eG4/inc/constructProtonBeamDump.hh"

#include <iostream>

#include "G4Color.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/ProtonBeamDump.hh"

#include "G4Helper/inc/VolumeInfo.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"

namespace mu2e {
  void constructProtonBeamDump(const VolumeInfo& parent, const SimpleConfig& config) {

    const ProtonBeamDump& det = *(GeomHandle<ProtonBeamDump>());    

    MaterialFinder materialFinder(config);
    
    VolumeInfo logicalEnclosure = nestBox("ProtonBeamDump",
					  det.enclosureHalfSize(), 
					  materialFinder.get("protonBeamDump.material.air"),
					  &det.enclosureRotationInMu2e(), // assume the parent is not rotated
					  det.enclosureCenterInMu2e() - parent.centerInMu2e(),
					  parent, 0, config.getBool("protonBeamDump.logicalEnclosureVisible"),
					  G4Colour::Grey(),
					  false, true, true, true
					  );

  }
}
