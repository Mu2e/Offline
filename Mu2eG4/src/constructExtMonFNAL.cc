// Andrei Gaponenko, 2011

#include "Mu2eG4/inc/constructExtMonFNAL.hh"

#include <iostream>

#include "G4Color.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

#include "G4Helper/inc/VolumeInfo.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"

#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"

namespace mu2e {
  void constructExtMonFNAL(const VolumeInfo& parent,
                           const CLHEP::HepRotation &parentRotationInMu2e,
                           const CLHEP::HepRotation *rotationInParent,
                           const SimpleConfig& config) {

    bool const forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    GeomHandle<ExtMonFNAL::ExtMon> det;
    GeomHandle<ProtonBeamDump> dump;

    MaterialFinder materialFinder(config);
    G4VSensitiveDetector* emSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::ExtMonFNAL());

    const VolumeInfo room = nestBox("ExtMonFNALRoom",
                                    det->roomHalfSize(),
                                    materialFinder.get("extmon_fnal.roomMaterialName"),
                                    rotationInParent,
                                    parentRotationInMu2e*(det->roomCenterInMu2e() - parent.centerInMu2e()),
                                    parent,
                                    0,
                                    config.getBool("extmon_fnal.roomVisible"),
                                    G4Colour::Grey(),
                                    config.getBool("extmon_fnal.roomSolid"),
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );

    const VolumeInfo detector = nestBox("ExtMonFNALDetector",
                                        det->detectorHalfSize(),
                                        materialFinder.get("extmon_fnal.roomMaterialName"),
                                        &det->detectorRotationInRoom(),
                                        det->detectorCenterInRoom(),
                                        room, 0,
                                        config.getBool("extmon_fnal.detectorBoxVisible"),
                                        G4Colour::Red(),
                                        config.getBool("extmon_fnal.detectorBoxSolid"),
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
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
                                  config.getBool("extmon_fnal.detectorPlaneVisible"),
                                  G4Colour::Magenta(),
                                  config.getBool("extmon_fnal.detectorPlaneSolid"),
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck
                                  );

      vplane.logical->SetSensitiveDetector(emSD);
    }
  }
}
