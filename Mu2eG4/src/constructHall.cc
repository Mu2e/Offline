//
// Free function to create the hall walls and hall interior inside the earthen overburden.
//
//
// Original author KLG based on Mu2eWorld constructHall
//
// Notes:
// Construct the earthen overburden

// Mu2e includes
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeneralUtilities/inc/OrientationResolver.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/NotchManager.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eHallGeom/inc/Mu2eHall.hh"
#include "Mu2eG4/inc/constructHall.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/nestBox.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4ExtrudedSolid.hh"
#include "G4GenericTrap.hh"
#include "G4RotationMatrix.hh"
#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4TwoVector.hh"
#include "CLHEP/Vector/Rotation.h"
#include "G4NistManager.hh"

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
    NotchManager notchMgr;
    notchMgr.loadNotches(config);

    constructSolids( config, hallInfo, building->getBldgSolids(), horizontalConcreteRotation, notchMgr );
    constructSolids( config, hallInfo, building->getDirtSolids(), horizontalConcreteRotation, notchMgr );
    constructTrapSolids( config, hallInfo, building->getDirtTrapSolids(), horizontalConcreteRotation, notchMgr );

    return hallInfo;

  }

  //================================================================================
  void constructSolids( const SimpleConfig& config,
                        const VolumeInfo& hallInfo, 
			const std::map<std::string,ExtrudedSolid>& solidMap,
			const CLHEP::HepRotation& rot,
			const NotchManager& notchMgr) {
    
    //-----------------------------------------------------------------
    // Building and dirt volumes are extruded solids.
    //-----------------------------------------------------------------
    
    const auto& geoOptions         = art::ServiceHandle<GeometryService>()->geomOptions();
    const bool doSurfaceCheck      = geoOptions->doSurfaceCheck("HallAir"); 
    const bool forceAuxEdgeVisible = geoOptions->forceAuxEdgeVisible("HallAir"); 
    const bool placePV             = geoOptions->placePV("HallAir"); 

    OrientationResolver* OR = new OrientationResolver();

    // Loop over all volumes in the map
    for ( const auto& keyVolumePair : solidMap ) {

      const auto& volume = keyVolumePair.second;
      const auto& volName = keyVolumePair.first;
      
      geoOptions->loadEntry( config, volume.getName(), volume.getName() );

      if ( notchMgr.hasNotches( volName ) ) {
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
      
	// Now loop over the notches and subtract them from above
	// First, create the eventual solid
	G4SubtractionSolid* aSolid = 0;
	// Get the vector of notches
	vector<Notch> volNotches = notchMgr.getNotchVector(volName);
	for ( unsigned int iNotch = 0; iNotch < volNotches.size(); iNotch++ ) {
	  ostringstream notchName;
	  notchName << "Notch" << iNotch+1;
	  Notch tmpNotch = volNotches[iNotch];
	  vector<double> halfDims = tmpNotch.getDims();
	  G4Box* notchBox = new G4Box( notchName.str(), 
				       halfDims[0],halfDims[1],halfDims[2]);

	  CLHEP::HepRotation* notchRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	  OR->getRotationFromOrientation( *notchRotat, tmpNotch.getOrient());
				       
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

      } else {  // Now for walls without notches

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

      } // end else of if for has notches

    } // end loop over parts

    // clean up a bit
    if ( 0 != OR ) {
      delete OR;
      OR = 0;
    }
  } // end function def for constructSolids

  //================================================================================
  void constructTrapSolids( const SimpleConfig& config,
			    const VolumeInfo& hallInfo, 
			    const std::map<std::string,GenericTrap>& solidMap,
			    const CLHEP::HepRotation& rot,
			    const NotchManager& notchMgr) {
    
    //-----------------------------------------------------------------
    // Building and dirt volumes are generic trapezoids.
    //-----------------------------------------------------------------
    
    const auto& geoOptions         = art::ServiceHandle<GeometryService>()->geomOptions();
    const bool doSurfaceCheck      = geoOptions->doSurfaceCheck("HallAir"); 
    const bool forceAuxEdgeVisible = geoOptions->forceAuxEdgeVisible("HallAir"); 
    const bool placePV             = geoOptions->placePV("HallAir"); 

    OrientationResolver* OR = new OrientationResolver();

    // Loop over all volumes in the map
    for ( const auto& keyVolumePair : solidMap ) {

      const auto& volume = keyVolumePair.second;
      const auto& volName = keyVolumePair.first;
      
      geoOptions->loadEntry( config, volume.getName(), volume.getName() );
      const CLHEP::HepRotation vRot = rot*volume.getRotation();
      const G4RotationMatrix* vRotG4 = new G4RotationMatrix(vRot);

      if ( notchMgr.hasNotches( volName ) ) {
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
	vector<Notch> volNotches = notchMgr.getNotchVector(volName);
	for ( unsigned int iNotch = 0; iNotch < volNotches.size(); iNotch++ ) {
	  ostringstream notchName;
	  notchName << "Notch" << iNotch+1;
	  Notch tmpNotch = volNotches[iNotch];
	  vector<double> halfDims = tmpNotch.getDims();
	  G4Box* notchBox = new G4Box( notchName.str(), 
				       halfDims[0],halfDims[1],halfDims[2]);

	  CLHEP::HepRotation* notchRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	  OR->getRotationFromOrientation( *notchRotat, tmpNotch.getOrient());
				       
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

      } else {  // Now for volumes without notches

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

    // clean up a bit
    if ( 0 != OR ) {
      delete OR;
      OR = 0;
    }
  } // end function def for constructTrapSolids

} // end namespace mu2e
