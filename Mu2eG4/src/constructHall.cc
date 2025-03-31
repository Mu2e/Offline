//
// Free function to create the hall walls and hall interior inside the earthen overburden.
//
//
// Original author KLG based on Mu2eWorld constructHall
//
// Notes:
// Construct the earthen overburden

// Mu2e includes
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/GeneralUtilities/inc/OrientationResolver.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/NotchHoleManager.hh"
#include "Offline/GeometryService/inc/WorldG4.hh"
#include "Offline/Mu2eHallGeom/inc/Mu2eHall.hh"
#include "Offline/Mu2eG4/inc/constructHall.hh"
#include "Offline/Mu2eG4/inc/MaterialFinder.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"

// G4 includes
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4ExtrudedSolid.hh"
#include "Geant4/G4GenericTrap.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4Orb.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4TwoVector.hh"
#include "CLHEP/Vector/Rotation.h"
#include "Geant4/G4NistManager.hh"

// C++ includes
#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>

using namespace std;

namespace mu2e {

  VolumeInfo constructHall(const VolumeInfo& worldInfo, const SimpleConfig& config ) {

    MaterialFinder materialFinder( config );

    GeomHandle<WorldG4> world;
    GeomHandle<Mu2eHall> building;

    const auto& geoOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geoOptions->loadEntry( config, "hallFormalBox", "hall.formalBox" );
    geoOptions->loadEntry( config, "HallAir", "hall" );
    const bool isHallFormalBoxVisible = geoOptions->isVisible("hallFormalBox");
    const bool isHallFormalBoxSolid   = geoOptions->isSolid("hallFormalBox");
    const bool doSurfaceCheck         = geoOptions->doSurfaceCheck("HallAir");
    const bool forceAuxEdgeVisible    = geoOptions->forceAuxEdgeVisible("HallAir");
    const bool placePV                = geoOptions->placePV("HallAir");

    // The formal hall volume
    VolumeInfo hallInfo = nestBox( "HallAir",
                                   world->hallFormalHalfSize(),
                                   materialFinder.get("hall.insideMaterialName"),
                                   0,
                                   world->hallFormalCenterInWorld(),
                                   worldInfo,
                                   0,
                                   isHallFormalBoxVisible,
                                   G4Colour::Red(),
                                   isHallFormalBoxSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   );

    // Rotation is static because rotations are not copied into G4.
    static CLHEP::HepRotation horizontalConcreteRotation(CLHEP::HepRotation::IDENTITY);
    horizontalConcreteRotation.rotateX( 90*CLHEP::degree);
    horizontalConcreteRotation.rotateZ( 90*CLHEP::degree);

    // Allow notches/holes in building walls
    NotchHoleManager notchMgr;
    notchMgr.loadNotches(config);
    notchMgr.loadHoles(config);

    constructSolids( config, hallInfo, building->getBldgSolids(), horizontalConcreteRotation, notchMgr );
    constructSolids( config, hallInfo, building->getDirtSolids(), horizontalConcreteRotation, notchMgr );
    constructRotSolids( config, hallInfo, building->getRotSolids(), horizontalConcreteRotation, notchMgr );
    constructTrapSolids( config, hallInfo, building->getDirtTrapSolids(), horizontalConcreteRotation, notchMgr );

    return hallInfo;

  }

  //================================================================================
  void constructSolids( const SimpleConfig& config,
                        const VolumeInfo& hallInfo,
                        const std::map<std::string,ExtrudedSolid>& solidMap,
                        const CLHEP::HepRotation& rot,
                        const NotchHoleManager& notchMgr) {

    //-----------------------------------------------------------------
    // Building and dirt volumes are extruded solids.
    //-----------------------------------------------------------------

    const auto& geoOptions         = art::ServiceHandle<GeometryService>()->geomOptions();
    const bool doSurfaceCheck      = geoOptions->doSurfaceCheck("HallAir");
    const bool forceAuxEdgeVisible = geoOptions->forceAuxEdgeVisible("HallAir");
    const bool placePV             = geoOptions->placePV("HallAir");

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    OrientationResolver OR;

    // Loop over all volumes in the map
    for ( const auto& keyVolumePair : solidMap ) {

      const auto& volume = keyVolumePair.second;
      const auto& volName = keyVolumePair.first;

      geoOptions->loadEntry( config, volume.getName(), volume.getName() );

      if ( notchMgr.hasNotches( volName )||notchMgr.hasHoles( volName ) ) {
       const vector<Notch>& volNotches = notchMgr.getNotchVector(volName);
       const vector<Hole>& volHoles = notchMgr.getHoleVector(volName);

       // First do volumes with notches
       // Make the VolumeInfo, without solid info
        VolumeInfo tmpVol(volume.getName(),
                          volume.getOffsetFromMu2eOrigin() - hallInfo.centerInMu2e(),
                          hallInfo.centerInWorld);
       // Make the main extruded solid from which notches will be subtracted
        G4ExtrudedSolid* aVol = new G4ExtrudedSolid(tmpVol.name,
                                                    volume.getVertices(),
                                                    volume.getYhalfThickness(),
                                                    G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);
       //-------------------------------------------------------------------------------------
       // Now loop over the notches and subtract them from above
       // First, create the eventual solid
       G4SubtractionSolid* aSolid = 0;
       //-------------------------------------------------------------------------------------
       // Get the vector of notches
        for ( unsigned int iNotch = 0; iNotch < volNotches.size(); iNotch++ ) {
          ostringstream notchName;
          notchName << "Notch" << iNotch+1;
          Notch tmpNotch = volNotches[iNotch];
          vector<double> halfDims = tmpNotch.getDims();
          G4Box* notchBox = new G4Box( notchName.str(),
                                       halfDims[0],halfDims[1],halfDims[2]);

          CLHEP::HepRotation* notchRotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
          OR.getRotationFromOrientation( *notchRotat, tmpNotch.getOrient());

          if ( 0 == aSolid ) {
            aSolid = new G4SubtractionSolid( tmpVol.name,
                                             aVol,
                                             notchBox,
                                             notchRotat,
                                             tmpNotch.getCenter() );
          } else {
            G4SubtractionSolid * bSolid = new G4SubtractionSolid
              ( tmpVol.name,
                aSolid,
                notchBox,
                notchRotat,
                tmpNotch.getCenter() );
            aSolid = bSolid;
          } // end if...else for first or later notch
        } // end loop over all notches
       //-----------------------------------------------------------------------------------
       //-------------------------------------------------------------------------------------
       // Get the vector of holes
        for ( unsigned int iHole = 0; iHole < volHoles.size(); iHole++ ) {
          ostringstream holeName;
          holeName << "hole" << iHole+1;
          Hole tmpHole = volHoles[iHole];
          double radius = tmpHole.getRad();
          double halfLength = tmpHole.getHalfLen();
          G4Tubs*holeTub = new G4Tubs(holeName.str(),
                                      0.0, //inner radius
                                      radius,//outer radius
                                      halfLength,  // Half Lengths are kept at least 1 mm longer.
                                      0.0, //angle_0
                                      CLHEP::twopi); //angle_span;

          CLHEP::HepRotation* holeRotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
          OR.getRotationFromOrientation( *holeRotat, tmpHole.getOrient());
          if ( 0 == aSolid ) {
            aSolid = new G4SubtractionSolid( tmpVol.name,
                                             aVol,
                                             holeTub,
                                             holeRotat,
                                             tmpHole.getCenter() );
          } else {
            G4SubtractionSolid * bSolid = new G4SubtractionSolid
              ( tmpVol.name,
                aSolid,
                holeTub,
                holeRotat,
                tmpHole.getCenter() );
            aSolid = bSolid;
          } // end if...else for first or later hole
        } // end loop over all holes
        //-----------------------------------------------------------------------------------

        tmpVol.solid = aSolid;

        finishNesting(tmpVol,
                      findMaterialOrThrow( volume.getMaterial() ),
                      &rot,
                      tmpVol.centerInParent,
                      hallInfo.logical,
                      0,
                      geoOptions->isVisible( volume.getName() ),
                      G4Colour::Grey(),
                      geoOptions->isSolid( volume.getName() ),
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck
                      );

      } else {  // Now for walls without notches or holes

        VolumeInfo tmpVol(volume.getName(),
                          volume.getOffsetFromMu2eOrigin() - hallInfo.centerInMu2e(),
                          hallInfo.centerInWorld);

        tmpVol.solid = new G4ExtrudedSolid(tmpVol.name,
                                           volume.getVertices(),
                                           volume.getYhalfThickness(),
                                           G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);

        finishNesting(tmpVol,
                      findMaterialOrThrow( volume.getMaterial() ),
                      &rot,
                      tmpVol.centerInParent,
                      hallInfo.logical,
                      0,
                      geoOptions->isVisible( volume.getName() ),
                      G4Colour::Grey(),
                      geoOptions->isSolid( volume.getName() ),
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck
                      );

      } // end else of if for has notches and holes

    } // end loop over parts

  } // end function def for constructSolids

  //================================================================================
  void constructRotSolids( const SimpleConfig& config,
                           const VolumeInfo& hallInfo,
                           const std::map<std::string,RotExtrudedSolid>& solidMap,
                           const CLHEP::HepRotation& rot,
                           const NotchHoleManager& notchMgr) {

    //-----------------------------------------------------------------
    // Some building and dirt volumes (mainly stairs) require to rotate extruded solids.
    //-----------------------------------------------------------------
    const auto& geoOptions         = art::ServiceHandle<GeometryService>()->geomOptions();
    const bool doSurfaceCheck      = geoOptions->doSurfaceCheck("HallAir");
    const bool forceAuxEdgeVisible = geoOptions->forceAuxEdgeVisible("HallAir");
    const bool placePV             = geoOptions->placePV("HallAir");

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    OrientationResolver OR;

    // Loop over all volumes in the map
    for ( const auto& keyVolumePair : solidMap ) {

      const auto& volume = keyVolumePair.second;
      const auto& volName = keyVolumePair.first;

      geoOptions->loadEntry( config, volume.getName(), volume.getName() );
      const CLHEP::HepRotation vRot = rot*volume.getRotation();
      const G4RotationMatrix* vRotG4 = reg.add(G4RotationMatrix(vRot));

      if ( notchMgr.hasNotches( volName )||notchMgr.hasHoles( volName ) ) {
       const vector<Notch>& volNotches = notchMgr.getNotchVector(volName);
       const vector<Hole>& volHoles = notchMgr.getHoleVector(volName);

       // First do volumes with notches
       // Make the VolumeInfo, without solid info
        VolumeInfo tmpVol(volume.getName(),
                          volume.getOffsetFromMu2eOrigin() - hallInfo.centerInMu2e(),
                          hallInfo.centerInWorld);
       // Make the main extruded solid from which notches will be subtracted
        G4ExtrudedSolid* aVol = new G4ExtrudedSolid(tmpVol.name,
                                                    volume.getVertices(),
                                                    volume.getYhalfThickness(),
                                                    G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);
       //-------------------------------------------------------------------------------------
       // Now loop over the notches and subtract them from above
       // First, create the eventual solid
       G4SubtractionSolid* aSolid = 0;
       //-------------------------------------------------------------------------------------
       // Get the vector of notches
        for ( unsigned int iNotch = 0; iNotch < volNotches.size(); iNotch++ ) {
          ostringstream notchName;
          notchName << "Notch" << iNotch+1;
          Notch tmpNotch = volNotches[iNotch];
          vector<double> halfDims = tmpNotch.getDims();
          G4Box* notchBox = new G4Box( notchName.str(),
                                       halfDims[0],halfDims[1],halfDims[2]);

          CLHEP::HepRotation* notchRotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
          OR.getRotationFromOrientation( *notchRotat, tmpNotch.getOrient());

          if ( 0 == aSolid ) {
            aSolid = new G4SubtractionSolid( tmpVol.name,
                                             aVol,
                                             notchBox,
                                             notchRotat,
                                             tmpNotch.getCenter() );
          } else {
            G4SubtractionSolid * bSolid = new G4SubtractionSolid
              ( tmpVol.name,
                aSolid,
                notchBox,
                notchRotat,
                tmpNotch.getCenter() );
            aSolid = bSolid;
          } // end if...else for first or later notch
        } // end loop over all notches
       //-----------------------------------------------------------------------------------
       //-------------------------------------------------------------------------------------
       // Get the vector of holes
        for ( unsigned int iHole = 0; iHole < volHoles.size(); iHole++ ) {
          ostringstream holeName;
          holeName << "hole" << iHole+1;
          Hole tmpHole = volHoles[iHole];
          double radius = tmpHole.getRad();
          double halfLength = tmpHole.getHalfLen();
          G4Tubs*holeTub = new G4Tubs(holeName.str(),
                                      0.0, //inner radius
                                      radius,//outer radius
                                      halfLength,  // Half Lengths are kept at least 1 mm longer.
                                      0.0, //angle_0
                                      CLHEP::twopi); //angle_span;

          CLHEP::HepRotation* holeRotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
          OR.getRotationFromOrientation( *holeRotat, tmpHole.getOrient());
          if ( 0 == aSolid ) {
            aSolid = new G4SubtractionSolid( tmpVol.name,
                                             aVol,
                                             holeTub,
                                             holeRotat,
                                             tmpHole.getCenter() );
          } else {
            G4SubtractionSolid * bSolid = new G4SubtractionSolid
              ( tmpVol.name,
                aSolid,
                holeTub,
                holeRotat,
                tmpHole.getCenter() );
            aSolid = bSolid;
          } // end if...else for first or later hole
        } // end loop over all holes
        //-----------------------------------------------------------------------------------

        tmpVol.solid = aSolid;

        finishNesting(tmpVol,
                      findMaterialOrThrow( volume.getMaterial() ),
                      vRotG4,
                      tmpVol.centerInParent,
                      hallInfo.logical,
                      0,
                      geoOptions->isVisible( volume.getName() ),
                      G4Colour::Grey(),
                      geoOptions->isSolid( volume.getName() ),
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck
                      );

      } else {  // Now for walls without notches or holes
        VolumeInfo tmpVol(volume.getName(),
                          volume.getOffsetFromMu2eOrigin() - hallInfo.centerInMu2e(),
                          hallInfo.centerInWorld);
        tmpVol.solid = new G4ExtrudedSolid(tmpVol.name,
                                           volume.getVertices(),
                                           volume.getYhalfThickness(),
                                           G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);

        finishNesting(tmpVol,
                      findMaterialOrThrow( volume.getMaterial() ),
                      vRotG4,
                      tmpVol.centerInParent,
                      hallInfo.logical,
                      0,
                      geoOptions->isVisible( volume.getName() ),
                      G4Colour::Grey(),
                      geoOptions->isSolid( volume.getName() ),
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck
                      );
      } // end else of if for has notches and holes

    } // end loop over parts

  } // end function def for constructSolids

  //================================================================================
  void constructTrapSolids( const SimpleConfig& config,
                            const VolumeInfo& hallInfo,
                            const std::map<std::string,GenericTrap>& solidMap,
                            const CLHEP::HepRotation& rot,
                            const NotchHoleManager& notchMgr) {

    //-----------------------------------------------------------------
    // Building and dirt volumes are generic trapezoids.
    //-----------------------------------------------------------------

    const auto& geoOptions         = art::ServiceHandle<GeometryService>()->geomOptions();
    const bool doSurfaceCheck      = geoOptions->doSurfaceCheck("HallAir");
    const bool forceAuxEdgeVisible = geoOptions->forceAuxEdgeVisible("HallAir");
    const bool placePV             = geoOptions->placePV("HallAir");

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    OrientationResolver OR;

    // Loop over all volumes in the map
    for ( const auto& keyVolumePair : solidMap ) {

      const auto& volume = keyVolumePair.second;
      const auto& volName = keyVolumePair.first;

      geoOptions->loadEntry( config, volume.getName(), volume.getName() );
      const CLHEP::HepRotation vRot = rot*volume.getRotation();
      const G4RotationMatrix* vRotG4 = reg.add(G4RotationMatrix(vRot));

      if ( notchMgr.hasNotches( volName )||notchMgr.hasHoles( volName ) ) {
        // First do volumes with notches

        // Make the VolumeInfo, without solid info
        VolumeInfo tmpVol(volume.getName(),
                          volume.getOffsetFromMu2eOrigin() - hallInfo.centerInMu2e(),
                          hallInfo.centerInWorld);

        // Make the main extruded solid from which notches will be subtracted
        G4GenericTrap* aVol = new G4GenericTrap(tmpVol.name,
                                                volume.getYhalfThickness(),
                                                volume.getVertices()
                                                );

        // Now loop over the notches and subtract them from above
        // First, create the eventual solid
        G4SubtractionSolid* aSolid = 0;
        // Get the vector of notches
        const vector<Notch>& volNotches = notchMgr.getNotchVector(volName);
        for ( unsigned int iNotch = 0; iNotch < volNotches.size(); iNotch++ ) {
          ostringstream notchName;
          notchName << "Notch" << iNotch+1;
          Notch tmpNotch = volNotches[iNotch];
          vector<double> halfDims = tmpNotch.getDims();
          G4Box* notchBox = new G4Box( notchName.str(),
                                       halfDims[0],halfDims[1],halfDims[2]);

          CLHEP::HepRotation* notchRotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
          OR.getRotationFromOrientation( *notchRotat, tmpNotch.getOrient());

          if ( 0 == aSolid ) {
            aSolid = new G4SubtractionSolid( tmpVol.name,
                                             aVol,
                                             notchBox,
                                             notchRotat,
                                             tmpNotch.getCenter() );
          } else {
            G4SubtractionSolid * bSolid = new G4SubtractionSolid
              ( tmpVol.name,
                aSolid,
                notchBox,
                notchRotat,
                tmpNotch.getCenter() );
            aSolid = bSolid;
          } // end if...else for first or later notch
        } // end loop over all notches

        //-------------------------holes---------------------------------------------
        // Get the vector of holes
        const vector<Hole>& volHoles = notchMgr.getHoleVector(volName);
        for ( unsigned int iHole = 0; iHole < volHoles.size(); iHole++ ) {
          ostringstream holeName;
          holeName << "Hole" << iHole+1;
          Hole tmpHole = volHoles[iHole];
          double radius = tmpHole.getRad();
          double halfLength = tmpHole.getHalfLen();
          G4Tubs*holeTub = new G4Tubs(holeName.str(),
                                      0.0, //inner radius
                                      radius,//outer radius
                                      halfLength, //height
                                      0.0, //angle_0
                                      CLHEP::twopi); //angle_span

          CLHEP::HepRotation* holeRotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
          OR.getRotationFromOrientation( *holeRotat, tmpHole.getOrient());

          if ( 0 == aSolid ) {
            aSolid = new G4SubtractionSolid( tmpVol.name,
                                             aVol,
                                             holeTub,
                                             holeRotat,
                                             tmpHole.getCenter() );
          } else {
            G4SubtractionSolid * bSolid = new G4SubtractionSolid
              ( tmpVol.name,
                aSolid,
                holeTub,
                holeRotat,
                tmpHole.getCenter() );
            aSolid = bSolid;
          } // end if...else for first or later hole
        } // end loop over all holes

        tmpVol.solid = aSolid;

        finishNesting(tmpVol,
                      findMaterialOrThrow( volume.getMaterial() ),
                      vRotG4,
                      tmpVol.centerInParent,
                      hallInfo.logical,
                      0,
                      geoOptions->isVisible( volume.getName() ),
                      G4Colour::Grey(),
                      geoOptions->isSolid( volume.getName() ),
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck
                      );

      } else {  // Now for volumes without notches or holes

        VolumeInfo tmpVol(volume.getName(),
                          volume.getOffsetFromMu2eOrigin() - hallInfo.centerInMu2e(),
                          hallInfo.centerInWorld);

        tmpVol.solid = new G4GenericTrap(tmpVol.name,
                                         volume.getYhalfThickness(),
                                         volume.getVertices()
                                         );

        auto vertices = volume.getVertices();
        finishNesting(tmpVol,
                      findMaterialOrThrow( volume.getMaterial() ),
                      // 0,
                      vRotG4,
                      tmpVol.centerInParent,
                      hallInfo.logical,
                      0,
                      geoOptions->isVisible( volume.getName() ),
                      G4Colour::Grey(),
                      geoOptions->isSolid( volume.getName() ),
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck
                      );

      } // end else of if for has notches

    } // end loop over parts

   } // end function def for constructTrapSolids

} // end namespace mu2e
