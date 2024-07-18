// Sam Fine, 2024

#include "Offline/Mu2eG4/inc/constructExtMonFNAL.hh"

#include <vector>
#include <numbers>
#include <math.h>

#include "Geant4/G4Color.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4Orb.hh"
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4ExtrudedSolid.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "cetlib_except/exception.h"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"

#include "Offline/Mu2eHallGeom/inc/Mu2eHall.hh"

#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/Mu2eG4Helper/inc/AntiLeakRegistry.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeomPrimitives/inc/TubsParams.hh"
#include "Offline/GeomPrimitives/inc/ExtrudedSolid.hh"
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
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "extMonFNAL", "extMonFNAL" );

    const bool forceAuxEdgeVisible  = geomOptions->forceAuxEdgeVisible("extMonFNAL");
    const bool doSurfaceCheck       = geomOptions->doSurfaceCheck("extMonFNAL");
    const bool placePV              = geomOptions->placePV("extMonFNAL");

    //----------------------------------------------------------
    // Pixel Chiller Mother Volume

    G4RotationMatrix* pixelChillerRotationInParentG4 = reg.add(G4RotationMatrix(pixelChillerRotationInParent));
    std::vector<double> extMonFNALPixelChillerHalfSize;
    config.getVectorDouble("extMonFNAL.pixelChiller.halfSize", extMonFNALPixelChillerHalfSize);
    const std::string pixelChillerNum = "pixelChiller"+std::to_string(chillerNum);

    VolumeInfo pixelChillerMother =
      nestBox(pixelChillerNum,
              extMonFNALPixelChillerHalfSize,
              findMaterialOrThrow("G4_AIR"),
              pixelChillerRotationInParentG4,
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

    std::vector<double> pixelChillerTopPCBParams = { extMonFNALPixelChillerHalfSize[0] - 2 * (pixelChillerPCBThickness + PCBDistanceFromWall),
                                                     pixelChillerPCBThickness,
                                                     extMonFNALPixelChillerHalfSize[2] - 2 * (pixelChillerPCBThickness + PCBDistanceFromWall)
                                                   };

    double pixelChillerTopBotPCBOffset = extMonFNALPixelChillerHalfSize[1] - PCBDistanceFromWall - 0.5 * pixelChillerPCBThickness;

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
                                                           extMonFNALPixelChillerHalfSize[1] - 2 * (pixelChillerPCBThickness + PCBDistanceFromWall),
                                                           extMonFNALPixelChillerHalfSize[2] - 2 * (pixelChillerPCBThickness + PCBDistanceFromWall)
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

    std::vector<double> pixelChillerEastSidePCBParams = { extMonFNALPixelChillerHalfSize[0] - 2 * (pixelChillerPCBThickness + PCBDistanceFromWall),
                                                          extMonFNALPixelChillerHalfSize[1] - 2 * (pixelChillerPCBThickness + PCBDistanceFromWall),
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

//================================================================
void constructElectronicsRackExtMonFNAL( const VolumeInfo& parent,
                                         const CLHEP::Hep3Vector& eRackCenterInParent,
                                         const CLHEP::HepRotation& eRackRotationInParent,
                                         const SimpleConfig& config
                                       )
  {
    MaterialFinder materialFinder(config);
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "extMonFNAL", "extMonFNAL" );

    const bool forceAuxEdgeVisible  = geomOptions->forceAuxEdgeVisible("extMonFNAL");
    const bool doSurfaceCheck       = geomOptions->doSurfaceCheck("extMonFNAL");
    const bool placePV              = geomOptions->placePV("extMonFNAL");

    //----------------------------------------------------------
    // Electronics Rack Mother Volume

    G4RotationMatrix* eRackRotationInParentG4 = reg.add(G4RotationMatrix(eRackRotationInParent));
    std::vector<double> extMonFNALElectronicsRackHalfSize;
    config.getVectorDouble("extMonFNAL.electronicsRack.halfSize", extMonFNALElectronicsRackHalfSize);

      const VolumeInfo electronicsRackMother =
        nestBox( "ExtMonFNALElectronicsRackMother",
                 extMonFNALElectronicsRackHalfSize,
                 findMaterialOrThrow("G4_AIR"),
                 eRackRotationInParentG4,
                 eRackCenterInParent,
                 parent, 0,
                 geomOptions->isVisible("ExtMonFNALElectronicsRack"),
                 G4Colour::Magenta(),
                 geomOptions->isSolid("ExtMonFNALElectronicsRack"),
                 forceAuxEdgeVisible,
                 placePV,
                 doSurfaceCheck
               );

    //----------------------------------------------------------
    // Electronics Rack Steel Case

    VolumeInfo electronicsRackSteelCase( "ExtMonFNALElectronicsRackSteelCase",
                                         CLHEP::Hep3Vector(0,0,0),
                                         electronicsRackMother.centerInWorld
                                       );

    G4Box* electronicsRackSteelCaseBox = new G4Box ( "electronicsRackSteelCaseBox",
                                                     extMonFNALElectronicsRackHalfSize[0],
                                                     extMonFNALElectronicsRackHalfSize[1],
                                                     extMonFNALElectronicsRackHalfSize[2]
                                                   );

    double electronicsRackCaseThickness = config.getDouble("extMonFNAL.electronicsRackCaseThickness");

    G4Box* electronicsRackCaseSubtractionBox = new G4Box ( "electronicsRackSteelCaseSubtractionBox",
                                                            extMonFNALElectronicsRackHalfSize[0] - electronicsRackCaseThickness,
                                                            extMonFNALElectronicsRackHalfSize[1] - electronicsRackCaseThickness,
                                                            extMonFNALElectronicsRackHalfSize[2] - electronicsRackCaseThickness
                                                          );

    electronicsRackSteelCase.solid = new G4SubtractionSolid( "electronicsRackSteelCase",
                                                             electronicsRackSteelCaseBox,
                                                             electronicsRackCaseSubtractionBox
                                                           );

    finishNesting(electronicsRackSteelCase,
                  materialFinder.get("extMonFNAL.electronicsRack.caseMaterialName"),
                  0,
                  CLHEP::Hep3Vector(0,0,0),
                  electronicsRackMother.logical,
                  0,
                  geomOptions->isVisible( "electronicsRackSteelCase"),
                  G4Colour::Red() ,
                  geomOptions->isSolid( "electronicsRackSteelCase"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    //----------------------------------------------------------
    // Electronics Rack Circuit Boards

    double electronicsRackPCBThickness = config.getDouble("extMonFNAL.electronicsRack.PCBThickness");
    double PCBDistanceFromWall = config.getDouble("extMonFNAL.electronicsRack.PCBDistanceFromWall");

    std::vector<double> electronicsRackTopPCBParams = { extMonFNALElectronicsRackHalfSize[0] - 2 * (electronicsRackPCBThickness + PCBDistanceFromWall),
                                                        electronicsRackPCBThickness,
                                                        extMonFNALElectronicsRackHalfSize[2] - 2 * (electronicsRackPCBThickness + PCBDistanceFromWall)
                                                      };

    double electronicsRackTopBotPCBOffset = extMonFNALElectronicsRackHalfSize[1] - PCBDistanceFromWall - 0.5 * electronicsRackPCBThickness;

    CLHEP::Hep3Vector electronicsRackTopPCBCenter ( 0, electronicsRackTopBotPCBOffset, 0 );

    VolumeInfo electronicsRackPCBTop =
      nestBox("electronicsRackTopPCB",
              electronicsRackTopPCBParams,
              materialFinder.get("extMonFNAL.electronicsRack.PCBMaterialName"),
              0,
              electronicsRackTopPCBCenter,
              electronicsRackMother.logical,
              0,
              geomOptions->isVisible("electronicsRack"),
              G4Colour::Red(),
              geomOptions->isSolid("electronicsRack"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    CLHEP::Hep3Vector electronicsRackBottomPCBCenter ( 0, -electronicsRackTopBotPCBOffset, 0 );

    VolumeInfo electronicsRackPCBBottom =
      nestBox("electronicsRackBottomPCB",
              electronicsRackTopPCBParams,
              materialFinder.get("extMonFNAL.electronicsRack.PCBMaterialName"),
              0,
              electronicsRackBottomPCBCenter,
              electronicsRackMother.logical,
              0,
              geomOptions->isVisible("electronicsRack"),
              G4Colour::Red(),
              geomOptions->isSolid("electronicsRack"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    std::vector<double> electronicsRackNorthSidePCBParams = { electronicsRackPCBThickness,
                                                              extMonFNALElectronicsRackHalfSize[1] - 2 * (electronicsRackPCBThickness + PCBDistanceFromWall),
                                                              extMonFNALElectronicsRackHalfSize[2] - 2 * (electronicsRackPCBThickness + PCBDistanceFromWall)
                                                            };

    double electronicsRackNSPCBOffset = extMonFNALElectronicsRackHalfSize[0] - PCBDistanceFromWall - 0.5 * electronicsRackPCBThickness;

    CLHEP::Hep3Vector electronicsRackNorthPCBCenter ( electronicsRackNSPCBOffset, 0, 0 );

    VolumeInfo electronicsRackNorthSidePCB =
      nestBox("electronicsRackNorthSidePCB",
              electronicsRackNorthSidePCBParams,
              materialFinder.get("extMonFNAL.electronicsRack.PCBMaterialName"),
              0,
              electronicsRackNorthPCBCenter,
              electronicsRackMother.logical,
              0,
              geomOptions->isVisible("electronicsRack"),
              G4Colour::Red(),
              geomOptions->isSolid("electronicsRack"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    CLHEP::Hep3Vector electronicsRackSouthPCBCenter ( - electronicsRackNSPCBOffset, 0, 0 );

    VolumeInfo electronicsRackSouthSidePCB =
      nestBox("electronicsRackSouthSidePCB",
              electronicsRackNorthSidePCBParams,
              materialFinder.get("extMonFNAL.electronicsRack.PCBMaterialName"),
              0,
              electronicsRackSouthPCBCenter,
              electronicsRackMother.logical,
              0,
              geomOptions->isVisible("electronicsRack"),
              G4Colour::Red(),
              geomOptions->isSolid("electronicsRack"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    std::vector<double> electronicsRackEastSidePCBParams = { extMonFNALElectronicsRackHalfSize[0] - 2 * (electronicsRackPCBThickness + PCBDistanceFromWall),
                                                             extMonFNALElectronicsRackHalfSize[1] - 2 * (electronicsRackPCBThickness + PCBDistanceFromWall),
                                                             electronicsRackPCBThickness
                                                           };

    double electronicsRackEWPCBOffset = extMonFNALElectronicsRackHalfSize[2] - PCBDistanceFromWall - 0.5 * electronicsRackPCBThickness;

    CLHEP::Hep3Vector electronicsRackEastPCBCenter ( 0, 0, electronicsRackEWPCBOffset );

    VolumeInfo electronicsRackEastSidePCB =
      nestBox("electronicsRackEastSidePCB",
              electronicsRackEastSidePCBParams,
              materialFinder.get("extMonFNAL.electronicsRack.PCBMaterialName"),
              0,
              electronicsRackEastPCBCenter,
              electronicsRackMother.logical,
              0,
              geomOptions->isVisible("electronicsRack"),
              G4Colour::Red(),
              geomOptions->isSolid("electronicsRack"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    CLHEP::Hep3Vector electronicsRackWestPCBCenter ( 0, 0, -electronicsRackEWPCBOffset );

    VolumeInfo electronicsRackWestSidePCB =
      nestBox("electronicsRackWestSidePCB",
              electronicsRackEastSidePCBParams,
              materialFinder.get("extMonFNAL.electronicsRack.PCBMaterialName"),
              0,
              electronicsRackWestPCBCenter,
              electronicsRackMother,
              0,
              geomOptions->isVisible("electronicsRack"),
              G4Colour::Red(),
              geomOptions->isSolid("electronicsRack"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    //----------------------------------------------------------
    // Electronics Rack Stacked Ciruit Boards

    double numLayers = config.getDouble("extMonFNAL.electronicsRack.numLayers");
    std::vector<double> electronicsRackStackPCBParams = { extMonFNALElectronicsRackHalfSize[0] - 2 * (electronicsRackPCBThickness + PCBDistanceFromWall),
                                                          electronicsRackPCBThickness,
                                                          extMonFNALElectronicsRackHalfSize[2] - 2 * (electronicsRackPCBThickness + PCBDistanceFromWall)
                                                        };

    double electronicsRackStackSize = extMonFNALElectronicsRackHalfSize[1] - PCBDistanceFromWall - 0.5 * electronicsRackPCBThickness;

    for (int i = 1; i <= numLayers; i++) {

    CLHEP::Hep3Vector electronicsRackStackPCBCenter ( 0, -electronicsRackStackSize + 2 * i * electronicsRackStackSize / (numLayers + 1), 0 );

    VolumeInfo electronicsRackPCBStack =
      nestBox("electronicsRackStack"+std::to_string(i)+"PCB",
              electronicsRackStackPCBParams,
              materialFinder.get("extMonFNAL.electronicsRack.PCBMaterialName"),
              0,
              electronicsRackStackPCBCenter,
              electronicsRackMother.logical,
              0,
              geomOptions->isVisible("electronicsRack"),
              G4Colour::Red(),
              geomOptions->isSolid("electronicsRack"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );
    }
  }

  //================================================================
  void constructExtMonFNALInfrastructure(const VolumeInfo& pixelChillerParent,
                                         const CLHEP::HepRotation& pixelChillerParentRotationInMu2e,
                                         const VolumeInfo& mainParent,
                                         const CLHEP::HepRotation& mainParentRotationInMu2e,
                                         const SimpleConfig& config)
  {
    MaterialFinder materialFinder(config);
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "extMonFNAL", "extMonFNAL" );
    GeomHandle<ExtMonFNALBuilding> emfb;
    GeomHandle<Mu2eHall> hall;

    //----------------------------------------------------------
    // Construct Pixel Chillers in Detector Room

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
        chillerCenter += CLHEP::Hep3Vector(emfb->detectorRoomHalfSize()[0] - extMonFNALPixelChillerHalfSize[0]  - distanceToSideWall,
                                           -emfb->detectorRoomHalfSize()[1] + extMonFNALPixelChillerHalfSize[1] + distanceToFloor,
                                           emfb->detectorRoomHalfSize()[2] - extMonFNALPixelChillerHalfSize[2] - distanceToTranverseWall
                                          );
      } else if (cornerReference == "SE") {
        chillerCenter += CLHEP::Hep3Vector( -emfb->detectorRoomHalfSize()[0] + extMonFNALPixelChillerHalfSize[0] + distanceToSideWall,
                                            -emfb->detectorRoomHalfSize()[1] + extMonFNALPixelChillerHalfSize[1] + distanceToFloor,
                                            emfb->detectorRoomHalfSize()[2] - extMonFNALPixelChillerHalfSize[2] - distanceToTranverseWall
                                          );
      } else if (cornerReference == "SW") {
        chillerCenter += CLHEP::Hep3Vector( -emfb->detectorRoomHalfSize()[0] + extMonFNALPixelChillerHalfSize[0] + distanceToSideWall,
                                            -emfb->detectorRoomHalfSize()[1] + extMonFNALPixelChillerHalfSize[1] + distanceToFloor,
                                            -emfb->detectorRoomHalfSize()[2] + extMonFNALPixelChillerHalfSize[2] + distanceToTranverseWall
                                          );
      } else if (cornerReference == "NW") {
        chillerCenter += CLHEP::Hep3Vector( emfb->detectorRoomHalfSize()[0] - extMonFNALPixelChillerHalfSize[0]  - distanceToSideWall,
                                            -emfb->detectorRoomHalfSize()[1] + extMonFNALPixelChillerHalfSize[1] + distanceToFloor,
                                            -emfb->detectorRoomHalfSize()[2] + extMonFNALPixelChillerHalfSize[2] + distanceToTranverseWall
                                          );
      } else {
        throw cet::exception("CONFIG") << "Error: constructExtMonFNALInfrastructure() cannot parse cornerReference = " << cornerReference<<"\n";
      }

      constructPixelChillerExtMonFNAL(i,
                                      pixelChillerParent,
                                      chillerCenter,
                                      CLHEP::HepRotation::IDENTITY,
                                      config);
    }

    //----------------------------------------------------------
    // Construct Electronics Rack

    double eRackDistanceToFloor = config.getDouble("extMonFNAL.room.electronicsRack.distanceToFloor");
    double eRackDistanceToSideWall = config.getDouble("extMonFNAL.room.electronicsRack.distanceToSideWall");
    double eRackDistanceToTranverseWall = config.getDouble("extMonFNAL.room.electronicsRack.distanceToTranverseWall");
    std::vector<double> extMonFNALElectronicsRackHalfSize;
    config.getVectorDouble("extMonFNAL.electronicsRack.halfSize", extMonFNALElectronicsRackHalfSize);

    CLHEP::Hep3Vector eRackCenter ( -emfb->detectorRoomHalfSize()[0] + extMonFNALElectronicsRackHalfSize[0] + eRackDistanceToSideWall,
                                    -emfb->detectorRoomHalfSize()[1] + extMonFNALElectronicsRackHalfSize[1] + eRackDistanceToFloor,
                                    -emfb->detectorRoomHalfSize()[2] + extMonFNALElectronicsRackHalfSize[2] + eRackDistanceToTranverseWall);

    constructElectronicsRackExtMonFNAL( pixelChillerParent,
                                        eRackCenter,
                                        CLHEP::HepRotation::IDENTITY,
                                        config);

    //----------------------------------------------------------
    // Construct Pixel Chiller Along ExtMon Room Wall

    const std::string wallName = config.getString("extMonFNAL.hall.chiller1.wallName");
    std::vector<int> chillerRefVertices;
    config.getVectorInt("extMonFNAL.hall.chiller1.refVertices", chillerRefVertices);
    std::vector<double> extMonFNALPixelChillerHalfSize;
    config.getVectorDouble("extMonFNAL.pixelChiller.halfSize", extMonFNALPixelChillerHalfSize);
    const double distanceToFloor = config.getDouble("extMonFNAL.hall.chiller1.distanceToFloor");
    const double distanceToWall  = config.getDouble("extMonFNAL.hall.chiller1.distanceToWall");
    const double offsetAlongWall = config.getDouble("extMonFNAL.hall.chiller1.offsetAlongWall");
    const ExtrudedSolid& extMonRoomWall = hall->getBldgSolid(wallName);
    const auto& extMonRoomWallVerticies = extMonRoomWall.getVertices();
    const double halfWallHeight = extMonRoomWall.getYhalfThickness();

    using CLHEP::Hep2Vector;
    Hep2Vector corner1 (extMonRoomWallVerticies[chillerRefVertices[0]].y(), extMonRoomWallVerticies[chillerRefVertices[0]].x());
    Hep2Vector corner2 (extMonRoomWallVerticies[chillerRefVertices[1]].y(), extMonRoomWallVerticies[chillerRefVertices[1]].x());
    Hep2Vector wallVector = corner2 - corner1;
    Hep2Vector wallChillerXZ = corner1;
    wallChillerXZ -= Hep2Vector(mainParent.centerInMu2e().x(), mainParent.centerInMu2e().z());
    wallChillerXZ += Hep2Vector(extMonRoomWall.getOffsetFromMu2eOrigin().x(), extMonRoomWall.getOffsetFromMu2eOrigin().z());
    Hep2Vector wallVectorUnit = wallVector.unit();
    wallChillerXZ += offsetAlongWall * wallVectorUnit;
    wallChillerXZ -= distanceToWall * wallVectorUnit.orthogonal();
    wallChillerXZ -= extMonFNALPixelChillerHalfSize[2] * wallVectorUnit.orthogonal();
    wallChillerXZ += extMonFNALPixelChillerHalfSize[0] * wallVectorUnit;
    const double chillerYCoord = extMonRoomWall.getOffsetFromMu2eOrigin().y()
                                 - mainParent.centerInMu2e().y()
                                 - halfWallHeight + extMonFNALPixelChillerHalfSize[1] + distanceToFloor; //changes chiller's height off the floor

    double wallVectorX = wallVector.y();
    double wallVectorY = wallVector.x();
    double rotAngle = -M_PI * 0.5 - std::atan2 ( wallVectorY, wallVectorX );
    CLHEP::Hep3Vector chillerCoord ( wallChillerXZ.x(), chillerYCoord, wallChillerXZ.y());

    CLHEP::HepRotation *psideChillerRot = reg.add(new CLHEP::HepRotation());
    CLHEP::HepRotation& sideChillerRot = *psideChillerRot;
    psideChillerRot->rotateY( rotAngle );

    constructPixelChillerExtMonFNAL(5,
                                    mainParent,
                                    chillerCoord,
                                    sideChillerRot,
                                    config);

    //----------------------------------------------------------
    // Construct Pixel Chiller Along Coll2Shielding Wall

    Hep2Vector coll2WallCorner1 (emfb->coll2ShieldingOutline()[2].y(), emfb->coll2ShieldingOutline()[2].x());
    Hep2Vector coll2WallCorner2 (emfb->coll2ShieldingOutline()[3].y(), emfb->coll2ShieldingOutline()[3].x());
    Hep2Vector coll2WallVector = coll2WallCorner2 - coll2WallCorner1;
    Hep2Vector coll2WallChillerXZ = coll2WallCorner1;
    coll2WallChillerXZ -= Hep2Vector(mainParent.centerInMu2e().x(), mainParent.centerInMu2e().z());
    coll2WallChillerXZ += Hep2Vector(emfb->coll2ShieldingCenterInMu2e().x(), emfb->coll2ShieldingCenterInMu2e().z());
    Hep2Vector coll2WallVectorUnit = coll2WallVector.unit();
    coll2WallChillerXZ += offsetAlongWall * coll2WallVectorUnit;
    coll2WallChillerXZ -= distanceToWall * coll2WallVectorUnit.orthogonal();
    coll2WallChillerXZ -= extMonFNALPixelChillerHalfSize[2] * coll2WallVectorUnit.orthogonal();
    coll2WallChillerXZ += extMonFNALPixelChillerHalfSize[0] * coll2WallVectorUnit;
    const double coll2WallChillerYCoord = emfb->coll2ShieldingCenterInMu2e().y()
                                 - mainParent.centerInMu2e().y()
                                 - halfWallHeight + extMonFNALPixelChillerHalfSize[1] + distanceToFloor; //changes chiller's height off the floor

    CLHEP::Hep3Vector coll2WallChillerCoord ( coll2WallChillerXZ.x(), coll2WallChillerYCoord, coll2WallChillerXZ.y());
    double coll2WallVectorUnitX = coll2WallVectorUnit.y();
    double coll2WallVectorUnitY = coll2WallVectorUnit.x();
    double coll2WallRotAngle = -M_PI * 0.5 - std::atan2 ( coll2WallVectorUnitY, coll2WallVectorUnitX );

    CLHEP::HepRotation *psideColl2WallChillerRot = reg.add(new CLHEP::HepRotation());
    CLHEP::HepRotation& sideColl2WallChillerRot = *psideColl2WallChillerRot;
    psideColl2WallChillerRot->rotateY( coll2WallRotAngle );
    constructPixelChillerExtMonFNAL(6,
                                    mainParent,
                                    coll2WallChillerCoord,
                                    sideColl2WallChillerRot,
                                    config);
  } //constructExtMonFNALInfrastructure()
} //namespace mu2e
