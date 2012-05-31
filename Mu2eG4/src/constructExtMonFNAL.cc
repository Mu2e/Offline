//
// $Id: constructExtMonFNAL.cc,v 1.11 2012/05/31 17:09:13 genser Exp $
// $Author: genser $
// $Date: 2012/05/31 17:09:13 $
//
//
// Andrei Gaponenko, 2011

#include "Mu2eG4/inc/constructExtMonFNAL.hh"

#include <iostream>

#include "G4Color.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"
#include "G4ExtrudedSolid.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"

#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNALBuilding.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"

// FIXME: should not need WorldG4 here
#include "GeometryService/inc/WorldG4.hh"

namespace mu2e {
  void constructExtMonFNAL(const VolumeInfo& collimator1Parent,
                           const CLHEP::HepRotation& collimator1ParentRotationInMu2e,
                           const VolumeInfo& mainParent,
                           const CLHEP::HepRotation& mainParentRotationInMu2e,
                           const SimpleConfig& config)
  {
    const VolumeInfo roomAir = constructExtMonFNALBuilding(collimator1Parent,
                                                           collimator1ParentRotationInMu2e,
                                                           mainParent,
                                                           mainParentRotationInMu2e,
                                                           config);

    bool const forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
    bool const doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");
    bool const placePV             = true;

    MaterialFinder materialFinder(config);

    GeomHandle<ExtMonFNAL::ExtMon> extmon;
    GeomHandle<ExtMonFNALBuilding> emfb;

    G4VSensitiveDetector* emSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::ExtMonFNAL());


    //----------------------------------------------------------------
    // The detector enclosure

    // detectorRotationInRoom = roomRotationInMu2e.inverse() * detectorRotationInMu2e
    // We need the inverse: detRotInMu2e.inv() * roomRotationInMu2e
    static const CLHEP::HepRotation detectorRotationInRoomInv(
                                                              extmon->detectorRotationInMu2e().inverse()
                                                              * emfb->roomRotationInMu2e()
                                                              );

    const CLHEP::Hep3Vector detectorCenterInRoom(emfb->roomRotationInMu2e().inverse()*(extmon->detectorCenterInMu2e() - roomAir.centerInMu2e()));

    VolumeInfo detector = nestBox("ExtMonFNALDetector",
                                  extmon->detectorHalfSize(),
                                  materialFinder.get("extMonFNAL.room.materialName"),
                                  &detectorRotationInRoomInv,
                                  detectorCenterInRoom,
                                  roomAir, 0,
                                  config.getBool("extMonFNAL.detectorBoxVisible"),
                                  G4Colour::Red(),
                                  config.getBool("extMonFNAL.detectorBoxSolid"),
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck
                                  );

    // FIXME: we should not need to correct the wrong information
    detector.centerInWorld = GeomHandle<WorldG4>()->mu2eOriginInWorld() + extmon->detectorCenterInMu2e();

    //----------------------------------------------------------------
    // The sensor planes

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
                                  config.getBool("extMonFNAL.detectorPlaneVisible"),
                                  G4Colour::Magenta(),
                                  config.getBool("extMonFNAL.detectorPlaneSolid"),
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck
                                  );

      vplane.logical->SetSensitiveDetector(emSD);

    } // for()

    //----------------------------------------------------------------
    // The Virtual Detectors around ExtMonFNAL

    if(true) {
      const int verbosityLevel = config.getInt("vd.verbosityLevel");

      bool vdIsVisible         = config.getBool("vd.visible");
      bool vdIsSolid           = config.getBool("vd.solid");
      bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
      bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");
      bool const placePV       = true;

      int const nSurfaceCheckPoints = 100000; // for a more thorrow check due to the vd shape

      G4Material* vacuumMaterial     = materialFinder.get("toyDS.insideMaterialName");

      G4VSensitiveDetector* vdSD = G4SDManager::GetSDMpointer()->
        FindSensitiveDetector(SensitiveDetectorName::VirtualDetector());

      GeomHandle<VirtualDetector> vdg;

      // VDs around the detector
      for(int vdId = VirtualDetectorId::EMFDetectorEntrance; vdId <= VirtualDetectorId::EMFDetectorExit; ++vdId) {
        if( vdg->exist(vdId) ) {
          if ( verbosityLevel > 0) {
            std::cout<<__func__<<" constructing "<<VirtualDetector::volumeName(vdId)<<std::endl;
          }

          CLHEP::Hep3Vector centerInRoom = detectorCenterInRoom
            + detectorRotationInRoomInv.inverse()
            * CLHEP::Hep3Vector(0,
                                0,
                                (vdId == VirtualDetectorId::EMFDetectorEntrance ? 1 : -1)
                                * (extmon->detectorHalfSize()[2] + vdg->getHalfLength())
                                );


          std::vector<double> hlen(3);
          hlen[0] = config.getDouble("extMonFNAL.vd.halfdx");
          hlen[1] = config.getDouble("extMonFNAL.vd.halfdy");
          hlen[2] = vdg->getHalfLength();

          VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                      hlen,
                                      vacuumMaterial,
                                      &detectorRotationInRoomInv,
                                      centerInRoom,
                                      roomAir,
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
      } // for(vdId-1)

      //----------------------------------------------------------------
      for(int vdId = VirtualDetectorId::EMFC2Entrance; vdId <= VirtualDetectorId::EMFC2Exit; ++vdId) {
        if( vdg->exist(vdId) ) {
          if ( verbosityLevel > 0) {
            std::cout <<__func__<<" constructing "<<VirtualDetector::volumeName(vdId)<<std::endl;
          }

          std::vector<double> hlen(3);
          hlen[0] = emfb->coll2ShieldingHalfSize()[0];
          hlen[1] = emfb->coll2ShieldingHalfSize()[1];
          hlen[2] = vdg->getHalfLength();
          const CLHEP::Hep3Vector vdCenterInMu2e =
            emfb->coll2ShieldingCenterInMu2e()
            + emfb->coll2ShieldingRotationInMu2e()*CLHEP::Hep3Vector
            (0, 0,
             ((vdId == VirtualDetectorId::EMFC2Entrance) ? +1 : -1)*(emfb->coll2ShieldingHalfSize()[2] + hlen[2])
             );

          static const CLHEP::HepRotation vdRotInRoomInv
            (emfb->coll2ShieldingRotationInMu2e().inverse() * emfb->roomRotationInMu2e());

          VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                      hlen,
                                      vacuumMaterial,
                                      &vdRotInRoomInv,
                                      emfb->roomRotationInMu2e().inverse()*(vdCenterInMu2e - roomAir.centerInMu2e()),
                                      roomAir,
                                      vdId,
                                      vdIsVisible,
                                      G4Color::Red(),
                                      vdIsSolid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      false);

          // vd are very thin, a more thorough check is needed
          doSurfaceCheck && vdInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);

          vdInfo.logical->SetSensitiveDetector(vdSD);
        }
      } // for(vdId-2)

    } // VD block

  } // constructExtMonFNAL()
} // namespace mu2e
