// Sam Fine, 2024

#include "Offline/Mu2eG4/inc/constructExtMonFNAL.hh"

#include <vector>

#include "Geant4/G4Color.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4Box.hh"

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
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/Mu2eG4/inc/MaterialFinder.hh"

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
    // Pixel Chiller Volume

    std::vector<double> extMonFNALPixelChillerHalfSize;
    config.getVectorDouble("extMonFNAL.pixelChiller.halfSize", extMonFNALPixelChillerHalfSize);

    VolumeInfo pixelChiller =
      nestBox("pixelChiller"+std::to_string(chillerNum),
              extMonFNALPixelChillerHalfSize,
              materialFinder.get("extMonFNAL.pixelChiller.materialName"),
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

      // first coord corresponds to the northern axis
      // second cord corresponds to the 'z' axis (height off the floor of the room)
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
