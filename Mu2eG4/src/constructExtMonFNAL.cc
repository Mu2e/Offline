//
//
// Andrei Gaponenko, 2011

#include "Offline/Mu2eG4/inc/constructExtMonFNAL.hh"
#include "Offline/Mu2eG4/inc/constructExtMonFNALInfrastructure.hh"
#include "Offline/Mu2eG4/inc/constructExtMonFNALDetector.hh"

#include <iostream>

#include "Geant4/G4Color.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4ExtrudedSolid.hh"
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4Tubs.hh"
#include "Offline/Mu2eG4Helper/inc/AntiLeakRegistry.hh"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4/inc/nestExtrudedSolid.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/Mu2eG4/inc/MaterialFinder.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"

#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModuleIdConverter.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"
#include "Offline/Mu2eG4/inc/checkForOverlaps.hh"


//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {

  VolumeInfo constructExtMonFNALDetectorRoom(const VolumeInfo& parent,
                                             const CLHEP::HepRotation& parentRotationInMu2e,
                                             const SimpleConfig& config)
  {
    GeomHandle<ExtMonFNALBuilding> emfb;

    const std::string detectorRoomName = "ExtMonDetectorRoom";

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry(config, "extMonFNAL", "extMonFNAL");
    geomOptions->loadEntry(config, detectorRoomName, "extMonFNAL."+detectorRoomName);

    bool const isVisible = geomOptions->isVisible(detectorRoomName);
    bool const isSolid   = geomOptions->isSolid(detectorRoomName);
    bool const forceAuxEdgeVisible  = geomOptions->forceAuxEdgeVisible("extMonFNAL");
    bool const doSurfaceCheck       = geomOptions->doSurfaceCheck("extMonFNAL");
    bool const placePV              = geomOptions->placePV("extMonFNAL");

    //----------------------------------------------------------------
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    CLHEP::HepRotation *rotationInParentInv =
      reg.add(emfb->detectorRoomRotationInMu2e().inverse() * parentRotationInMu2e);

    const CLHEP::Hep3Vector refPointInParent(parentRotationInMu2e.inverse()*(emfb->detectorRoomCenterInMu2e() - parent.centerInMu2e()));

    //----------------------------------------------------------------

    static CLHEP::HepRotation collimator2ParentRotationInMu2e = emfb->coll2ShieldingRotationInMu2e();

    CLHEP::HepRotation *subCylinderRotation = reg.add(emfb->collimator2RotationInMu2e().inverse() * emfb->detectorRoomRotationInMu2e() );

    CLHEP::Hep3Vector subCylOffsetInParent = emfb->detectorRoomRotationInMu2e().inverse() * (emfb->collimator2CenterInMu2e() - emfb->detectorRoomCenterInMu2e());

    G4Tubs* subCylinder = new G4Tubs("detectorRoomSubtractionCylinder",
                                     0.*CLHEP::mm,
                                     emfb->collimator2().shotLinerOuterRadius(),
                                     0.5*emfb->collimator2().length(),
                                     0,
                                     CLHEP::twopi
                                     );

    VolumeInfo room(detectorRoomName,
                    refPointInParent,
                    parent.centerInWorld
                    );

    G4Box* roomBox = new G4Box("roomBox",
                               emfb->detectorRoomHalfSize()[0],
                               emfb->detectorRoomHalfSize()[1],
                               emfb->detectorRoomHalfSize()[2]
                               );

    room.solid = new G4SubtractionSolid(detectorRoomName,
                                        roomBox,
                                        subCylinder,
                                        subCylinderRotation,
                                        subCylOffsetInParent
                                        );

    finishNesting(room,
                  findMaterialOrThrow("G4_AIR"),
                  rotationInParentInv,
                  refPointInParent,
                  parent.logical,
                  0,
                  isVisible,
                  G4Colour::White(),
                  isSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    return room;
  }

  //================================================================
  void constructExtMonFNAL(const VolumeInfo& collimator1Parent,
                           const CLHEP::HepRotation& collimator1ParentRotationInMu2e,
                           const VolumeInfo& mainParent,
                           const CLHEP::HepRotation& mainParentRotationInMu2e,
                           const SimpleConfig& config)
  {
    constructExtMonFNALBuilding(collimator1Parent,
                                collimator1ParentRotationInMu2e,
                                mainParent,
                                mainParentRotationInMu2e,
                                config);

    VolumeInfo detectorRoom = constructExtMonFNALDetectorRoom(mainParent,
                                                              mainParentRotationInMu2e,
                                                              config);

    GeomHandle<ExtMonFNAL::ExtMon> extmon;
    GeomHandle<ExtMonFNALBuilding> emfb;

    constructExtMonFNALInfrastructure(detectorRoom,
                                      emfb->detectorRoomRotationInMu2e(),
                                      mainParent,
                                      mainParentRotationInMu2e,
                                      config);

    constructExtMonFNALDetector(detectorRoom,
                                emfb->detectorRoomRotationInMu2e(),
                                config);

  } // constructExtMonFNAL()
} // namespace mu2e
