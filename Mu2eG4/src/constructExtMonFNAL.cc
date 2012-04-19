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
#include "GeometryService/inc/VirtualDetector.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"
#include "Mu2eG4/inc/VirtualDetectorSD.hh"

namespace mu2e {
  void constructExtMonFNAL(const VolumeInfo& parent,
                           const CLHEP::HepRotation &parentRotationInMu2e,
                           const SimpleConfig& config) {

    bool const forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    GeomHandle<ExtMonFNAL::ExtMon> extmon;
    GeomHandle<ProtonBeamDump> dump;

    MaterialFinder materialFinder(config);
    G4VSensitiveDetector* emSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::ExtMonFNAL());


    // Room rotation in Mu2e == dump->enclosureRotationInMu2e()
    // So rotation in parent is parentRotationInMu2e.inverse() * dump->enclosureRotationInMu2e()
    //
    // finishNesting() uses the backwards interpretation of rotations, so it's
    // more convenient to store the inverse of that:
    static CLHEP::HepRotation roomRotationInParentInv(dump->enclosureRotationInMu2e().inverse()*parentRotationInMu2e);

    const VolumeInfo room = nestBox("ExtMonFNALRoom",
                                    extmon->roomHalfSize(),
                                    materialFinder.get("extmon_fnal.roomMaterialName"),
                                    &roomRotationInParentInv,
                                    parentRotationInMu2e.inverse()*(extmon->roomCenterInMu2e() - parent.centerInMu2e()),
                                    parent,
                                    0,
                                    config.getBool("extmon_fnal.roomVisible"),
                                    G4Colour::Grey(),
                                    config.getBool("extmon_fnal.roomSolid"),
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );


    // detectorRotationInRoom = roomRotationInMu2e.inverse() * detectorRotationInMu2e
    // We need the inverse: detRotInMu2e.inv() * roomRot
    static CLHEP::HepRotation detectorRotationInRoomInv(
                                                        extmon->detectorRotationInMu2e().inverse()
                                                        * dump->enclosureRotationInMu2e()
                                                        );

    const VolumeInfo detector = nestBox("ExtMonFNALDetector",
                                        extmon->detectorHalfSize(),
                                        materialFinder.get("extmon_fnal.roomMaterialName"),
                                        &detectorRotationInRoomInv, // det->detectorRotationInRoom(),
                                        extmon->detectorCenterInRoom(),
                                        room, 0,
                                        config.getBool("extmon_fnal.detectorBoxVisible"),
                                        G4Colour::Red(),
                                        config.getBool("extmon_fnal.detectorBoxSolid"),
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                        );

    for(unsigned iplane = 0; iplane < extmon->nplanes(); ++iplane) {
      std::ostringstream oss;
      oss<<"EMFSensor"<<iplane;

      VolumeInfo vplane = nestBox(oss.str(),
                                  extmon->sensorHalfSize(iplane),
                                  findMaterialOrThrow("G4_Si"),
                                  0,
                                  extmon->sensorOffsetInParent(iplane),
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

    } // for()




      //----------------------------------------------------------------
      // Construct the Virtual Detectors around ExtMonFNAL

    if(true) {
      int static const verbosityLevel = config.getInt("vd.verbosityLevel",0);

      bool vdIsVisible         = config.getBool("vd.visible",true);
      bool vdIsSolid           = config.getBool("vd.solid",true);
      bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
      bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
      bool const placePV       = true;

      int const nSurfaceCheckPoints = 100000; // for a more thorrow check due to the vd shape

      G4Material* vacuumMaterial     = materialFinder.get("toyDS.insideMaterialName");

      G4VSensitiveDetector* vdSD = G4SDManager::GetSDMpointer()->
        FindSensitiveDetector(SensitiveDetectorName::VirtualDetector());

      GeomHandle<VirtualDetector> vdg;

      for(int vdId = VirtualDetectorId::EMFDetectorEntrance; vdId <= VirtualDetectorId::EMFDetectorExit; ++vdId) {
        if( vdg->exist(vdId) ) {
          if ( verbosityLevel > 0) {
            std::cout<<__func__<<" constructing "<<VirtualDetector::volumeName(vdId)<<std::endl;
          }

          CLHEP::Hep3Vector centerInRoom = extmon->detectorCenterInRoom()
            + detectorRotationInRoomInv.inverse()
            * CLHEP::Hep3Vector(0,
                                0,
                                (vdId == VirtualDetectorId::EMFDetectorEntrance ? 1 : -1)
                                * (extmon->detectorHalfSize()[2] + vdg->getHalfLength())
                                );


          std::vector<double> hlen(3);
          hlen[0] = config.getDouble("extmon_fnal.vd.halfdx");
          hlen[1] = config.getDouble("extmon_fnal.vd.halfdy");
          hlen[2] = vdg->getHalfLength();

          VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                      hlen,
                                      vacuumMaterial,
                                      &detectorRotationInRoomInv,
                                      centerInRoom,
                                      room,
                                      vdId,
                                      vdIsVisible,
                                      G4Color::Red(),
                                      vdIsSolid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      false
                                      );

          // vd are very thin, a more thorough check is needed
          doSurfaceCheck && vdInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);

          vdInfo.logical->SetSensitiveDetector(vdSD);
        }
      }
    } // VD block

  } // constructExtMonFNAL()
} // namespace mu2e
