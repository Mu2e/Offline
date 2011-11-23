// Andrei Gaponenko, 2011

#include "Mu2eG4/inc/constructExtMonFNAL.hh"

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

#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"

namespace mu2e {
  void constructExtMonFNAL(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<ExtMonFNAL::ExtMon> det;
    GeomHandle<ProtonBeamDump> dump;

    MaterialFinder materialFinder(config);
    G4VSensitiveDetector* emSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::ExtMonFNAL());
    
    const VolumeInfo room = nestBox("ExtMonFNALRoom",
				    det->roomHalfSize(), 
				    materialFinder.get("extmon_fnal.roomMaterialName"),
				    
				    // rotated in the same way as the beam dump
				    &dump->enclosureRotationInMu2e(),
				    
				    det->roomCenterInMu2e() - parent.centerInMu2e(),
				    
				    parent, 0, config.getBool("extmon_fnal.roomVisible"),
				    G4Colour::Red(),
				    false, true, true, true
				    );
    
    const VolumeInfo detector = nestBox("ExtMonFNALDetector",
					det->detectorHalfSize(), 
					materialFinder.get("extmon_fnal.roomMaterialName"),
					&det->detectorRotationInRoom(),
					det->detectorCenterInRoom(),
					room, 0, false/*visible*/,
					G4Colour::Green(),
					false, true, true, true
					);
    
    for(unsigned iplane = 0; iplane < det->nplanes(); ++iplane) {
      std::ostringstream oss;
      oss<<"EMFSensor"<<iplane;

      VolumeInfo vplane = nestBox(oss.str(),
				  det->sensorHalfSize(iplane),
				  findMaterialOrThrow("G4_Si"),
				  0,
				  det->sensorOffsetInParent(iplane), 
				  detector,
				  0,
				  true, G4Colour::Gray(), true, true, true, true);
      
      vplane.logical->SetSensitiveDetector(emSD);
    }
  }
}
