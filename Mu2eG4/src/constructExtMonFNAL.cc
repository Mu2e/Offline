//
// $Id: constructExtMonFNAL.cc,v 1.20 2012/09/19 03:39:47 gandr Exp $
// $Author: gandr $
// $Date: 2012/09/19 03:39:47 $
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
#include "G4Helper/inc/G4Helper.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {

  //================================================================
  void constructExtMonFNALSensorStack(const ExtMonFNALSensor& sensor,
                                      const ExtMonFNALSensorStack& stack,
                                      const std::string& volNameSuffix,
                                      VirtualDetectorId::enum_type entranceVD,
                                      const VolumeInfo& parent,
                                      const CLHEP::HepRotation& parentRotationInMu2e,
                                      const SimpleConfig& config
                                      )
  {

    bool const forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
    bool const doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");
    bool const placePV             = true;

    MaterialFinder materialFinder(config);
    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

    G4VSensitiveDetector* emSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::ExtMonFNAL());

    //----------------------------------------------------------------
    // detectorRotationInRoom = roomRotationInMu2e.inverse() * detectorRotationInMu2e
    // We need the inverse: detRotInMu2e.inv() * roomRotationInMu2e

    CLHEP::HepRotation *stackRotationInRoomInv =
      reg.add(stack.rotationInMu2e().inverse() * parentRotationInMu2e);

    const CLHEP::HepRotation stackRotationInRoom(stackRotationInRoomInv->inverse());

    const CLHEP::Hep3Vector stackRefPointInRoom(parentRotationInMu2e.inverse()*(stack.refPointInMu2e() - parent.centerInMu2e()));

    //----------------------------------------------------------------
    // The sensor planes

    for(unsigned iplane = 0; iplane < stack.nplanes(); ++iplane) {
      std::ostringstream oss;
      oss<<"EMFSensor"<<volNameSuffix<<iplane;

      AGDEBUG("Constucting "<<oss.str()<<", sensor number "<<iplane + stack.planeNumberOffset());

      const CLHEP::Hep3Vector sensorCenterInRoom = stackRefPointInRoom + stackRotationInRoom*stack.sensorOffsetInStack(iplane);

      VolumeInfo vplane = nestBox(oss.str(),
                                  sensor.halfSize(),
                                  findMaterialOrThrow("G4_Si"),
                                  stackRotationInRoomInv,
                                  sensorCenterInRoom,
                                  parent,
                                  iplane + stack.planeNumberOffset(),
                                  config.getBool("extMonFNAL.sensorPlaneVisible"),
                                  G4Colour::Magenta(),
                                  config.getBool("extMonFNAL.sensorPlaneSolid"),
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck
                                  );

      vplane.logical->SetSensitiveDetector(emSD);

      // install passive readout material, adjacent to each sensor on
      // the downstream side

      std::ostringstream osr;
      osr<<"EMFReadout"<<volNameSuffix<<iplane;

      const CLHEP::Hep3Vector readoutCenterInRoom = stackRefPointInRoom
        + stackRotationInRoom*
        (stack.sensorOffsetInStack(iplane)
         + CLHEP::Hep3Vector(0,0, -(sensor.halfSize()[2]+stack.readout_halfdz()[iplane]))
         );

      std::vector<double> rhs(sensor.halfSize());
      rhs[2] = stack.readout_halfdz()[iplane];

      nestBox(osr.str(),
              rhs,
              findMaterialOrThrow("G4_Si"),
              stackRotationInRoomInv,
              readoutCenterInRoom,
              parent,
              0,
              config.getBool("extMonFNAL.readoutPlaneVisible"),
              G4Colour::Black(),
              config.getBool("extMonFNAL.readoutPlaneSolid"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );


    } // for()

    //----------------------------------------------------------------
    // Test material plates near the stack

    if(true) {
      AGDEBUG("constructing test materials: "<<stack.testMaterialNames().size()<<" plates");

      bool testMaterialVisible         = config.getBool("extMonFNAL.testMaterial.visible");
      bool testMaterialSolid           = config.getBool("extMonFNAL.testMaterial.solid");
      bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
      bool const placePV       = true;

      for(unsigned itest = 0; itest < stack.testMaterialNames().size(); ++itest) {

        CLHEP::Hep3Vector centerInRoom = stackRefPointInRoom
          + stackRotationInRoom
          * CLHEP::Hep3Vector(0,
                              0,
                              stack.distanceToTestMaterials()
                              + stack.testMaterialHalfSize()[2]
                              + itest * stack.testMaterialPitch()
                              );

        AGDEBUG("constructing test plate: "<<stack.testMaterialNames()[itest]);

        nestBox("ExtMonFNALTestMaterial_"+stack.testMaterialNames()[itest],
                stack.testMaterialHalfSize(),
                findMaterialOrThrow(stack.testMaterialNames()[itest]),
                stackRotationInRoomInv,
                centerInRoom,
                parent,
                0,
                testMaterialVisible,
                G4Color(1, 0.5, 1),
                testMaterialSolid,
                forceAuxEdgeVisible,
                placePV,
                false
                );

      } // for()
    } // test material block


    //----------------------------------------------------------------
    // VD block

    if(true) {

      const int verbosityLevel = config.getInt("vd.verbosityLevel");

      bool vdIsVisible         = config.getBool("vd.visible");
      bool vdIsSolid           = config.getBool("vd.solid");
      int const nSurfaceCheckPoints = 100000; // for a more thorrow check due to the vd shape

      MaterialFinder materialFinder(config);
      G4Material* vacuumMaterial     = materialFinder.get("toyDS.insideMaterialName");

      GeomHandle<VirtualDetector> vdg;

      G4VSensitiveDetector* vdSD = G4SDManager::GetSDMpointer()->
        FindSensitiveDetector(SensitiveDetectorName::VirtualDetector());

      for(int vdId = entranceVD; vdId <= 1 + entranceVD; ++vdId) {
        if( vdg->exist(vdId) ) {
          if ( verbosityLevel > 0) {
            std::cout<<__func__<<" constructing "<<VirtualDetector::volumeName(vdId)<<std::endl;
          }

          CLHEP::Hep3Vector centerInRoom = stackRefPointInRoom
            + stackRotationInRoom
            * CLHEP::Hep3Vector(0,
                                0,
                                (vdId == entranceVD ?
                                 (stack.sensorOffsetInStack(stack.nplanes()-1)[2]
                                  + sensor.halfSize()[2]
                                  + 2*stack.readout_halfdz()[stack.nplanes()-1]
                                  + vdg->getHalfLength()
                                  ) :
                                 (
                                  stack.sensorOffsetInStack(0)[2]
                                  - sensor.halfSize()[2]
                                  - 2*stack.readout_halfdz()[0]
                                  - vdg->getHalfLength()
                                  )
                                 )
                                );

          std::vector<double> hlen(3);
          hlen[0] = config.getDouble("extMonFNAL.vd.halfdx");
          hlen[1] = config.getDouble("extMonFNAL.vd.halfdy");
          hlen[2] = vdg->getHalfLength();

          VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                      hlen,
                                      vacuumMaterial,
                                      stackRotationInRoomInv,
                                      centerInRoom,
                                      parent,
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
    } // VD block

  }

  //================================================================
  void constructExtMonFNALVirtualDetectors(const VolumeInfo& roomAir,
                                           const CLHEP::HepRotation& parentRotationInMu2e,
                                           const SimpleConfig& config
                                           )
  {
    const int verbosityLevel = config.getInt("vd.verbosityLevel");

    bool vdIsVisible         = config.getBool("vd.visible");
    bool vdIsSolid           = config.getBool("vd.solid");
    bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
    bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");
    bool const placePV       = true;

    int const nSurfaceCheckPoints = 100000; // for a more thorrow check due to the vd shape

    MaterialFinder materialFinder(config);
    G4Material* vacuumMaterial     = materialFinder.get("toyDS.insideMaterialName");

    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

    G4VSensitiveDetector* vdSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::VirtualDetector());

    GeomHandle<VirtualDetector> vdg;
    GeomHandle<ExtMonFNALBuilding> emfb;

    //----------------------------------------------------------------
    const CLHEP::HepRotation* vdRotInRoomInv =
      reg.add(emfb->coll2ShieldingRotationInMu2e().inverse() * emfb->roomRotationInMu2e());

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

        VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                    hlen,
                                    vacuumMaterial,
                                    vdRotInRoomInv,
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
  }

  //================================================================
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

    GeomHandle<ExtMonFNAL::ExtMon> extmon;
    GeomHandle<ExtMonFNALBuilding> emfb;

    constructExtMonFNALSensorStack(extmon->sensor(),
                                   extmon->dn(),
                                   "Dn",
                                   VirtualDetectorId::EMFDetectorDnEntrance,
                                   roomAir,
                                   emfb->roomRotationInMu2e(),
                                   config);

    constructExtMonFNALSensorStack(extmon->sensor(),
                                   extmon->up(),
                                   "Up",
                                   VirtualDetectorId::EMFDetectorUpEntrance,
                                   roomAir,
                                   emfb->roomRotationInMu2e(),
                                   config);

    constructExtMonFNALMagnet(extmon->spectrometerMagnet(),
                              roomAir,
                              "spectrometer",
                              emfb->roomRotationInMu2e(),
                              config);


    constructExtMonFNALVirtualDetectors(roomAir, emfb->roomRotationInMu2e(), config);

  } // constructExtMonFNAL()
} // namespace mu2e
