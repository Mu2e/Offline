// Sam Fine, 2024

#include "Offline/Mu2eG4/inc/constructExtMonFNAL.hh"

#include <vector>
#include <numbers>

#include "Geant4/G4Color.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4SubtractionSolid.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "cetlib_except/exception.h"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"

#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/Mu2eG4Helper/inc/AntiLeakRegistry.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeomPrimitives/inc/TubsParams.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/Mu2eG4/inc/MaterialFinder.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"

namespace mu2e {

  //===============================================================
  void constructPixelChillerExtMonFNAL(const int chillerNum,
                                       const VolumeInfo& parent,
                                       const CLHEP::Hep3Vector& pixelChillerCenterInParent,
                                       const CLHEP::HepRotation& pixelChillerRotationInParent,
                                       const SimpleConfig& config
                                       )
  {
    MaterialFinder materialFinder(config);

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "extMonFNAL", "extMonFNAL" );

    const bool forceAuxEdgeVisible  = geomOptions->forceAuxEdgeVisible("extMonFNAL");
    const bool doSurfaceCheck       = geomOptions->doSurfaceCheck("extMonFNAL");
    const bool placePV              = geomOptions->placePV("extMonFNAL");

    //----------------------------------------------------------
    // Pixel Chiller Mother Volume

    std::vector<double> extMonFNALPixelChillerHalfSize;
    config.getVectorDouble("extMonFNAL.pixelChiller.halfSize", extMonFNALPixelChillerHalfSize);
    const std::string pixelChillerNum = "pixelChiller"+std::to_string(chillerNum);

    VolumeInfo pixelChillerMother =
      nestBox(pixelChillerNum,
              extMonFNALPixelChillerHalfSize,
              findMaterialOrThrow("G4_AIR"),
              0,
              pixelChillerCenterInParent,
              parent.logical,
              0,
              geomOptions->isVisible("pixelChiller"),
              G4Colour::Red(),
              geomOptions->isSolid("pixelChiller"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    //----------------------------------------------------------
    // Pixel Chiller Steel Case

    VolumeInfo pixelChillerSteelCase( "ExtMonFNAL"+pixelChillerNum+"SteelCase",
                                      CLHEP::Hep3Vector(0,0,0),
                                      pixelChillerMother.centerInWorld
                                    );

    G4Box* pixelChillerSteelCaseBox = new G4Box ( pixelChillerNum+"SteelCaseBox",
                                                  extMonFNALPixelChillerHalfSize[0],
                                                  extMonFNALPixelChillerHalfSize[1],
                                                  extMonFNALPixelChillerHalfSize[2]
                                                );

    double pixelChillerCaseThickness = config.getDouble("extMonFNAL.pixelChillerCaseThickness");

    G4Box* pixelChillerSteelCaseSubtractionBox = new G4Box ( pixelChillerNum+"SteelCaseSubtractionBox",
                                                             extMonFNALPixelChillerHalfSize[0] - pixelChillerCaseThickness,
                                                             extMonFNALPixelChillerHalfSize[1] - pixelChillerCaseThickness,
                                                             extMonFNALPixelChillerHalfSize[2] - pixelChillerCaseThickness
                                                           );

    pixelChillerSteelCase.solid = new G4SubtractionSolid( pixelChillerNum+"SteelCase",
                                                          pixelChillerSteelCaseBox,
                                                          pixelChillerSteelCaseSubtractionBox
                                                         );

    finishNesting(pixelChillerSteelCase,
                  materialFinder.get("extMonFNAL.pixelChiller.caseMaterialName"),
                  0,
                  CLHEP::Hep3Vector(0,0,0),
                  pixelChillerMother.logical,
                  0,
                  geomOptions->isVisible( "pixelChillerSteelCase"),
                  G4Colour::Red() ,
                  geomOptions->isSolid( "pixelChillerSteelCase"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    //----------------------------------------------------------
    // Pixel Chiller Circuit Boards

    double pixelChillerPCBThickness = config.getDouble("extMonFNAL.pixelChiller.PCBThickness");
    double PCBDistanceFromWall = config.getDouble("extMonFNAL.pixelChiller.PCBDistanceFromWall");

    std::vector<double> pixelChillerTopPCBParams = { extMonFNALPixelChillerHalfSize[0] - 2*(pixelChillerPCBThickness + PCBDistanceFromWall),
                                                     pixelChillerPCBThickness,
                                                     extMonFNALPixelChillerHalfSize[2] - 2*(pixelChillerPCBThickness + PCBDistanceFromWall)
                                                   };

    double pixelChillerTopBotPCBOffset = extMonFNALPixelChillerHalfSize[1] - PCBDistanceFromWall - 0.5* pixelChillerPCBThickness;

    CLHEP::Hep3Vector pixelChillerTopPCBCenter ( 0, pixelChillerTopBotPCBOffset, 0 );

    VolumeInfo pixelChillerPCBTop =
      nestBox(pixelChillerNum+"TopPCB",
              pixelChillerTopPCBParams,
              materialFinder.get("extMonFNAL.pixelChiller.PCBMaterialName"),
              0,
              pixelChillerTopPCBCenter,
              pixelChillerMother.logical,
              0,
              geomOptions->isVisible("pixelChiller"),
              G4Colour::Red(),
              geomOptions->isSolid("pixelChiller"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    CLHEP::Hep3Vector pixelChillerBottomPCBCenter ( 0, -pixelChillerTopBotPCBOffset, 0 );

    VolumeInfo pixelChillerPCBBottom =
      nestBox(pixelChillerNum+"BottomPCB",
              pixelChillerTopPCBParams,
              materialFinder.get("extMonFNAL.pixelChiller.PCBMaterialName"),
              0,
              pixelChillerBottomPCBCenter,
              pixelChillerMother.logical,
              0,
              geomOptions->isVisible("pixelChiller"),
              G4Colour::Red(),
              geomOptions->isSolid("pixelChiller"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    std::vector<double> pixelChillerNorthSidePCBParams = { pixelChillerPCBThickness,
                                                           extMonFNALPixelChillerHalfSize[1] - 2*(pixelChillerPCBThickness + PCBDistanceFromWall),
                                                           extMonFNALPixelChillerHalfSize[2] - 2*(pixelChillerPCBThickness + PCBDistanceFromWall)
                                                         };

    double pixelChillerNSPCBOffset = extMonFNALPixelChillerHalfSize[0] - PCBDistanceFromWall - 0.5 * pixelChillerPCBThickness;

    CLHEP::Hep3Vector pixelChillerNorthPCBCenter ( pixelChillerNSPCBOffset, 0, 0 );

    VolumeInfo pixelChillerNorthSidePCB =
      nestBox(pixelChillerNum+"NorthSidePCB",
              pixelChillerNorthSidePCBParams,
              materialFinder.get("extMonFNAL.pixelChiller.PCBMaterialName"),
              0,
              pixelChillerNorthPCBCenter,
              pixelChillerMother.logical,
              0,
              geomOptions->isVisible("pixelChiller"),
              G4Colour::Red(),
              geomOptions->isSolid("pixelChiller"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    CLHEP::Hep3Vector pixelChillerSouthPCBCenter ( - pixelChillerNSPCBOffset, 0, 0 );

    VolumeInfo pixelChillerSouthSidePCB =
      nestBox(pixelChillerNum+"SouthSidePCB",
              pixelChillerNorthSidePCBParams,
              materialFinder.get("extMonFNAL.pixelChiller.PCBMaterialName"),
              0,
              pixelChillerSouthPCBCenter,
              pixelChillerMother.logical,
              0,
              geomOptions->isVisible("pixelChiller"),
              G4Colour::Red(),
              geomOptions->isSolid("pixelChiller"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    std::vector<double> pixelChillerEastSidePCBParams = { extMonFNALPixelChillerHalfSize[0] - 2*(pixelChillerPCBThickness + PCBDistanceFromWall),
                                                          extMonFNALPixelChillerHalfSize[1] - 2*(pixelChillerPCBThickness + PCBDistanceFromWall),
                                                          pixelChillerPCBThickness
                                                        };

    double pixelChillerEWPCBOffset = extMonFNALPixelChillerHalfSize[2] - PCBDistanceFromWall - 0.5 * pixelChillerPCBThickness;

    CLHEP::Hep3Vector pixelChillerEastPCBCenter ( 0, 0, pixelChillerEWPCBOffset );

    VolumeInfo pixelChillerEastSidePCB =
      nestBox(pixelChillerNum+"EastSidePCB",
              pixelChillerEastSidePCBParams,
              materialFinder.get("extMonFNAL.pixelChiller.PCBMaterialName"),
              0,
              pixelChillerEastPCBCenter,
              pixelChillerMother.logical,
              0,
              geomOptions->isVisible("pixelChiller"),
              G4Colour::Red(),
              geomOptions->isSolid("pixelChiller"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    CLHEP::Hep3Vector pixelChillerWestPCBCenter ( 0, 0, -pixelChillerEWPCBOffset );

    VolumeInfo pixelChillerWestSidePCB =
      nestBox(pixelChillerNum+"WestSidePCB",
              pixelChillerEastSidePCBParams,
              materialFinder.get("extMonFNAL.pixelChiller.PCBMaterialName"),
              0,
              pixelChillerWestPCBCenter,
              pixelChillerMother,
              0,
              geomOptions->isVisible("pixelChiller"),
              G4Colour::Red(),
              geomOptions->isSolid("pixelChiller"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    //----------------------------------------------------------
    // Pixel Chiller Coolants

    double C5CoolantVolume = config.getDouble("extMonFNAL.pixelChiller.C5CoolantVolume");
    double freonVolume     = config.getDouble("extMonFNAL.pixelChiller.freonVolume");

    //Cylinder height is set as 4*radius, both are computed from volume
    double C5CoolantVolumeRadius = cbrt(C5CoolantVolume*CLHEP::liter / (4 * std::numbers::pi));
    double freonVolumeRadius = cbrt(freonVolume*CLHEP::liter / (4 * std::numbers::pi));

    CLHEP::Hep3Vector freonCenter ( freonVolumeRadius*2, 0, 0 );

    TubsParams pixelChillerFreonParams = { 0, freonVolumeRadius , 4*freonVolumeRadius };

    VolumeInfo pixelChillerFreon = nestTubs(pixelChillerNum+"Freon",
                                            pixelChillerFreonParams,
                                            materialFinder.get("extMonFNAL.pixelChiller.freonMaterial"),
                                            0,
                                            freonCenter,
                                            pixelChillerMother,
                                            0,
                                            geomOptions->isVisible("pixelChiller"),
                                            G4Colour::Red(),
                                            geomOptions->isSolid("pixelChiller"),
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                           );

    CLHEP::Hep3Vector CoolantVolCenter ( -C5CoolantVolumeRadius*2, 0, 0 );

    TubsParams pixelChillerC5CoolantParams = { 0, C5CoolantVolumeRadius , 4*C5CoolantVolumeRadius };

    VolumeInfo pixelChillerCoolant = nestTubs(pixelChillerNum+"Coolant",
                                              pixelChillerC5CoolantParams,
                                              materialFinder.get("extMonFNAL.pixelChiller.C5CoolantMaterial"),
                                              0,
                                              CoolantVolCenter,
                                              pixelChillerMother,
                                              0,
                                              geomOptions->isVisible("pixelChiller"),
                                              G4Colour::Red(),
                                              geomOptions->isSolid("pixelChiller"),
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                             );

  }

  void constructExtMonFNALInfrastructure(const VolumeInfo& pixelChillerParent,
                                         const CLHEP::HepRotation& pixelChillerParentRotationInMu2e,
                                         const VolumeInfo& mainParent,
                                         const CLHEP::HepRotation& mainParentRotationInMu2e,
                                         const SimpleConfig& config)
  {
    MaterialFinder materialFinder(config);
    GeomHandle<ExtMonFNALBuilding> emfb;

    int numChillers = config.getInt("extMonFNAL.numChillersInRoom");

    for (int i = 1; i <= numChillers; ++i) {

      std::string cornerReference = config.getString("extMonFNAL.room.chiller"+std::to_string(i)+".referenceCorner");
      double distanceToFloor = config.getDouble("extMonFNAL.room.chiller"+std::to_string(i)+".distanceToFloor");
      double distanceToSideWall = config.getDouble("extMonFNAL.room.chiller"+std::to_string(i)+".distanceToSideWall");
      double distanceToTranverseWall = config.getDouble("extMonFNAL.room.chiller"+std::to_string(i)+".distanceToTranverseWall");
      CLHEP::Hep3Vector chillerCenter (0,0,0);
      std::vector<double> extMonFNALPixelChillerHalfSize;
      config.getVectorDouble("extMonFNAL.pixelChiller.halfSize", extMonFNALPixelChillerHalfSize);

      // because of the ExtMon room volume rotation in Mu2e:
      // first coord corresponds to the northern axis
      // second coord corresponds to the vertical axis (height off the floor of the room)
      // third coord corresponds to the East/West axis

      if (cornerReference == "NE") {
        chillerCenter += CLHEP::Hep3Vector(emfb->detectorRoomHalfSize()[0] - extMonFNALPixelChillerHalfSize[0]*1.01  - distanceToSideWall,
                                           -emfb->detectorRoomHalfSize()[1] + extMonFNALPixelChillerHalfSize[1]*1.01 + distanceToFloor,
                                           emfb->detectorRoomHalfSize()[2] - extMonFNALPixelChillerHalfSize[2]*1.01 - distanceToTranverseWall
                                          );
      } else if (cornerReference == "SE") {
        chillerCenter += CLHEP::Hep3Vector( -emfb->detectorRoomHalfSize()[0] + extMonFNALPixelChillerHalfSize[0]*1.01 + distanceToSideWall,
                                            -emfb->detectorRoomHalfSize()[1] + extMonFNALPixelChillerHalfSize[1]*1.01 + distanceToFloor,
                                            emfb->detectorRoomHalfSize()[2] - extMonFNALPixelChillerHalfSize[2]*1.01 - distanceToTranverseWall
                                          );
      } else if (cornerReference == "SW") {
        chillerCenter += CLHEP::Hep3Vector( -emfb->detectorRoomHalfSize()[0] + extMonFNALPixelChillerHalfSize[0]*1.01 + distanceToSideWall,
                                            -emfb->detectorRoomHalfSize()[1] + extMonFNALPixelChillerHalfSize[1]*1.01 + distanceToFloor,
                                            -emfb->detectorRoomHalfSize()[2] + extMonFNALPixelChillerHalfSize[2]*1.01 + distanceToTranverseWall
                                          );
      } else if (cornerReference == "NW") {
        chillerCenter += CLHEP::Hep3Vector( emfb->detectorRoomHalfSize()[0] - extMonFNALPixelChillerHalfSize[0]*1.01  - distanceToSideWall,
                                            -emfb->detectorRoomHalfSize()[1] + extMonFNALPixelChillerHalfSize[1]*1.01 + distanceToFloor,
                                            -emfb->detectorRoomHalfSize()[2] + extMonFNALPixelChillerHalfSize[2]*1.01 + distanceToTranverseWall
                                          );
      } else {
        throw cet::exception("CONFIG")<< "Error: constructExtMonFNALInfrastructure() cannot parse cornerReference = "<<cornerReference<<"\n";
      }

      constructPixelChillerExtMonFNAL(i,
                                      pixelChillerParent,
                                      chillerCenter,
                                      CLHEP::HepRotation::IDENTITY,
                                      config);

    }

  } //constructExtMonFNALInfrastructure()
} //namespace mu2e
