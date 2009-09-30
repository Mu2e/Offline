//
// Construct the Mu2e G4 world and serve information about that world.
//
// $Id: Mu2eWorld.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//
//  Heirarchy is:
//  0      World (air)
//  1      Earthen Overburden
//  2      Concrete walls of the hall
//  3      Air inside the hall
//  4      Concrete shielding around the DS.
//  5      Air inside the concrete sheilding
//  6      Iron cosmic ray absorber.
//  7      Air inside the cosmic ray absorber
//  8      Effective volume representing the DS coils+cryostats.
//  9      Vacuum inside of the DS coils.
//
//  4      Effective volume representing the PS coils+cryostats.
//  5      Vacuum inside of the PS coils.
//
//  The Earth overburden is modelling in two parts: a box that extends
//  to the surface of the earth plus a cap above grade.  The cap is shaped
//  as a G4Paraboloid.
//

// C++ includes
#include <iostream>
#include <sstream>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/StrawPlacer.hh"
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"

// G4 includes
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Paraboloid.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"

using namespace std;

namespace mu2e {

  // Value used to request that all Mu2e specific materials be made.
  static const std::string DoAllValue = "DoAll";

  Mu2eWorld::Mu2eWorld()
    :  _cosmicReferencePoint(),
       _mu2eOrigin(),
       _info(),
       _detSolBField(),
       _usualRHS(),
       _exactHelix(),
       _chordFinder(),
       _fieldMgr(),
       _stepLimit(){
  }
  
  Mu2eWorld::~Mu2eWorld(){

    // Do not destruct the solids, logical volumes or physical volumes.
    // G4 looks after that itself.

  }

  // This is the callback used by G4.
  WorldInfo const* Mu2eWorld::construct(){
    
    // Get access to the master geometry system and its run time config.
    edm::Service<GeometryService> geom;
    SimpleConfig const& config = geom->config();
    _config = &config;
    
    // Construct a world with nothing in it.
    constructWorld(config);

    //constructTestWorld();
    return &_info;
  }

  void Mu2eWorld::constructWorld( SimpleConfig const& config ){

    // Dimensions and material of the world.
    vector<double> worldHLen;
    _config->getVectorDouble("world.halfLengths", worldHLen,3);
    setUnits ( worldHLen, mm );
    G4Material* worldMaterial = getMaterial("world.materialName");

    // A number of objects are referenced to the solenoids.
    double prodSolXoff = _config->getDouble("world.prodSolXoff") * mm;
    double detSolXoff  = -prodSolXoff;
    double dsz0        = _config->getDouble("toyDS.z0") * mm;

    // Half length of the detector solenoid.
    double dsHalfLength = _config->getDouble("toyDS.halfLength")* mm;
    
    // Construct the world volume.
    string worldName("World");
    _info.worldSolid = new G4Box( worldName, worldHLen[0], worldHLen[1], worldHLen[2] );
    _info.worldLog   = new G4LogicalVolume( _info.worldSolid, worldMaterial, worldName );
    _info.worldPhys  = new G4PVPlacement( 0, G4ThreeVector(), _info.worldLog, worldName, 
				     0, 0, 0);
    _info.worldLog->SetVisAttributes(G4VisAttributes::Invisible);

    // Get parameters related to the overall dimensions of the hall and to
    // the earthen overburden.
    double floorThick           = mm * _config->getDouble("hall.floorThick");
    double ceilingThick         = mm * _config->getDouble("hall.ceilingThick");
    double wallThick            = mm * _config->getDouble("hall.wallThick");
    double overburdenDepth      = mm * _config->getDouble("dirt.overburdenDepth");
    vector<double> hallInHLen;
    _config->getVectorDouble("hall.insideHalfLengths",hallInHLen,3);
    setUnits( hallInHLen, mm);

    // Derived parameters.
    G4Material* dirtMaterial = getMaterial("dirt.overburdenMaterialName");

    // Top of the floor in G4 world coordinates.
    double yFloor   = -worldHLen[1] + floorThick;

    // The height above the floor of the y origin of the Mu2e coordinate system.
    double yOriginHeight = _config->getDouble("world.mu2eOrigin.height" )*mm;

    // Position of the origin of the mu2e coordinate system
    _mu2eOrigin = G4ThreeVector( 
				_config->getDouble("world.mu2eOrigin.xoffset")*mm,
				yOriginHeight + yFloor,
				_config->getDouble("world.mu2eOrigin.zoffset")*mm
				);
    cout << "Mu2e Origin: " << _mu2eOrigin << endl;

    // Origin used to construct the MECO detector.
    _mu2eDetectorOrigin = _mu2eOrigin + G4ThreeVector( -3904., 0., 12000.);
    cout << "Mu2e Detector Origin: " << _mu2eDetectorOrigin << endl;

    // Bottom of the ceiling in G4 world coordinates.
    double yCeilingInSide = yFloor + 2.*hallInHLen[1];
    
    // Top of the ceiling in G4 world coordinates.
    double yCeilingOutside  = yCeilingInSide + ceilingThick;

    // Surface of the earth in G4 world coordinates.
    double ySurface  = yCeilingOutside + overburdenDepth;
    
    // Half length and y origin of the dirt box.
    double yLDirt = ( ySurface + worldHLen[1] )/2.;
    double y0Dirt = -worldHLen[1] + yLDirt;
    
    // Center of the dirt box, in the G4 world system.
    G4ThreeVector dirtOffset(0.,y0Dirt,0.);
    
    // Half lengths of the dirt box.
    double dirtHLen[3] = { worldHLen[0], yLDirt, worldHLen[2] };
    
    // Main body of dirt around the hall.
    VolumeInfo dirtInfo = nestBox( "DirtBody",
				   dirtHLen,
				   dirtMaterial,
				   0,
				   dirtOffset,
				   _info.worldLog,
				   0,
				   G4Colour::Magenta()
				   );
    
    // Dirt cap is modeled as a paraboloid.
    double capHalfHeight  = mm * _config->getDouble("dirt.capHalfHeight");
    double capBottomR     = mm * _config->getDouble("dirt.capBottomRadius");
    double capTopR        = mm * _config->getDouble("dirt.capTopRadius");

    // The top of the world.
    double yEverest = ySurface + 2.*capHalfHeight;

    // Selfconsistency check.
    if ( yEverest > 2.*worldHLen[1] ){
      throw cms::Exception("GEOM")
	<< "Top of the world is outside of the world volume! \n";
    }

    // Build the reference points that others will use.
    _cosmicReferencePoint = G4ThreeVector( 0., yEverest, 0.);
    cout << "Cosmic Ref: " << _cosmicReferencePoint << endl;

    // Construct the cap.
    string dirtCapName("DirtCap");
    G4Paraboloid* dirtCapSolid = new
      G4Paraboloid( dirtCapName, capHalfHeight, capTopR, capBottomR);
    
    G4LogicalVolume* dirtCapLog = new
      G4LogicalVolume( dirtCapSolid, 
		       dirtMaterial, 
		       dirtCapName);
    
    G4RotationMatrix* dirtCapRot = new G4RotationMatrix();
    dirtCapRot->rotateX( -90*degree);
    
    G4VPhysicalVolume* dirtCapPhys =  new
      G4PVPlacement( dirtCapRot, 
		     G4ThreeVector( detSolXoff, ySurface+capHalfHeight, dsz0+_mu2eOrigin.z()), 
		     dirtCapLog, 
		     dirtCapName, 
		     _info.worldLog, 
		     0, 
		     0);
    
    G4VisAttributes* dirtCapVisAtt = new G4VisAttributes(true, G4Colour::Green());
    dirtCapVisAtt->SetForceSolid(true);
    dirtCapLog->SetVisAttributes(dirtCapVisAtt);

    // Position of the center of the hall in the world volume.
    vector<double> hallPosition;
    _config->getVectorDouble("hall.offset", hallPosition,3);
    setUnits( hallPosition, mm);
    double hallY0 = yFloor + hallInHLen[1] + hallPosition[1];

    // Materials for the hall walls and the interior of the hall
    G4Material* wallMaterial = getMaterial("hall.wallMaterialName");
    G4Material* hallMaterial = getMaterial("hall.insideMaterialName");

    // Half lengths of the exterior of the concrete for the hall walls.
    double hallOutHLen[3] ={
      hallInHLen[0] + wallThick,
      hallInHLen[1] + ( ceilingThick + floorThick )/2.,
      hallInHLen[2] + wallThick
    };
    
    // Center of the concrete volume in the coordinate system of the dirt.
    G4ThreeVector wallOffset = 
      G4ThreeVector(hallPosition[0], hallY0, hallPosition[2]) - dirtOffset;

    // Origin of the hall air volume in the system of the hall concrete volume.
    G4ThreeVector hallOffset( 0., (floorThick-ceilingThick)/2., 0.);

    // Concrete walls of the hall.
    VolumeInfo wallInfo = nestBox( "HallWalls",
				   hallOutHLen,
				   wallMaterial,
				   0,
				   wallOffset,
				   dirtInfo.logical,
				   0,
				   G4Colour::Red()
				   );
    
    // Air volume inside of the hall.
    VolumeInfo hallInfo = nestBox( "HallAir",
				   hallInHLen,
				   hallMaterial,
				   0,
				   hallOffset,
				   wallInfo.logical,
				   0,
				   G4Colour::Red()
				   );
    
    // Concrete shield.
    double shieldConXSpace           = mm * _config->getDouble("shieldCon.xspace");
    double shieldConInsideHeight     = mm * _config->getDouble("shieldCon.insideHeight");
    double shieldConInsideHalfLength = mm * _config->getDouble("shieldCon.insideHalfLength");
    double shieldConThick            = mm * _config->getDouble("shieldCon.Thick");

    // The iron cosmic ray shield.
    double shieldFeThick     = mm * _config->getDouble("shieldFe.thick");
    double shieldFeOuterHW   = mm * _config->getDouble("shieldFe.outerHalfWidth");
    double shieldFeZExtra    = mm * _config->getDouble("shieldFe.zextra");
    
    // Materials for the above.
    G4Material* shieldConMaterial       = getMaterial("shieldCon.materialName");
    G4Material* shieldConInsideMaterial = getMaterial("shieldCon.insideMaterialName");
    G4Material* shieldFeMaterial        = getMaterial("shieldFe.materialName");
    G4Material* shieldFeInsideMaterial  = getMaterial("shieldFe.insideMaterialName");
    
    // Derived information for the concrete and Fe shields.
    double shieldConInsideHalfDim[3], shieldConOutsideHalfDim[3];
    double shieldFeInsideHalfDim[3], shieldFeOutsideHalfDim[3];
    
    shieldConOutsideHalfDim[0] = shieldFeOuterHW + shieldConXSpace + shieldConThick;
    shieldConOutsideHalfDim[1] = ( shieldConInsideHeight + shieldConThick)/2.;
    
    shieldConInsideHalfDim[0]  = shieldFeOuterHW + shieldConXSpace;
    shieldConInsideHalfDim[1]  = shieldConInsideHeight/2.;

    shieldFeOutsideHalfDim[0]  = shieldFeOuterHW;
    shieldFeOutsideHalfDim[1]  = shieldFeOuterHW;

    shieldFeInsideHalfDim[0]   = shieldFeOuterHW - shieldFeThick;
    shieldFeInsideHalfDim[1]   = shieldFeOuterHW - shieldFeThick;

     shieldFeInsideHalfDim[2]  = dsHalfLength + shieldFeZExtra;
    shieldFeOutsideHalfDim[2]  = shieldFeInsideHalfDim[2];

     shieldConInsideHalfDim[2] = shieldConInsideHalfLength;
    shieldConOutsideHalfDim[2] = shieldConInsideHalfLength + shieldConThick;

    // Position of concrete box inside the air volume of the hall.
    G4ThreeVector shieldConOffset = G4ThreeVector(
						  detSolXoff-hallPosition[0],
						  shieldConOutsideHalfDim[1] - hallInHLen[1],
						  dsz0+_mu2eOrigin.z()
						  );

    // Position of air inside concrete shield wrt concrete.
    double yoff2 = shieldConInsideHalfDim[1] - shieldConOutsideHalfDim[1];
    G4ThreeVector shieldConInsideOffset = G4ThreeVector(0.,yoff2,0.);

    // Position of the iron shield inside wrt the air inside the concrete.
    G4ThreeVector shieldFeOffset( 0., yOriginHeight-shieldConInsideHalfDim[1], 0. );

    // Concrete shield around DS.
    VolumeInfo shieldConInfo = nestBox( "ShieldConDS_01",
					shieldConOutsideHalfDim,
					shieldConMaterial,
					0,
					shieldConOffset,
					hallInfo.logical,
					0,
					G4Colour::Blue()
					);

    // Air between the concrete and Fe shields.
    VolumeInfo shieldConInsideInfo = nestBox( "ShieldConDS_01_AIR",
					      shieldConInsideHalfDim,
					      shieldConInsideMaterial,
					      0,
					      shieldConInsideOffset,
					      shieldConInfo.logical,
					      0,
					      G4Colour::Blue()
					      );
    
    // Fe shield around DS.
    VolumeInfo shieldFeInfo = nestBox( "ShieldFe_01",
				       shieldFeOutsideHalfDim,
				       shieldFeMaterial,
				       0,
				       shieldFeOffset,
				       shieldConInsideInfo.logical,
				       0,
				       G4Colour::Green()
				       );
    
    // Air between the Fe shield and the DS cryostat.
    VolumeInfo shieldFeInsideInfo = nestBox( "ShieldFe_AIR_01",
					     shieldFeInsideHalfDim,
					     shieldFeInsideMaterial,
					     0,
					     G4ThreeVector(),
					     shieldFeInfo.logical,
					     0,
					     G4Colour::Green()
					     );
    
    
    double detSolCoilParams[5] = { 
      _config->getDouble("toyDS.rOut"      ) * mm,
      _config->getDouble("toyDS.rIn"       ) * mm,
      _config->getDouble("toyDS.halfLength") * mm,
      0.,
      2.*M_PI
    };
    double detSolVacParams[5] = { 
      0. * mm,
      detSolCoilParams[0],
      _config->getDouble("toyDS.halfLengthVac") * mm,
      0.,
      2.*M_PI
    };
    G4Material* detSolCoilMaterial = getMaterial("toyDS.materialName");
    G4Material* detSolVacMaterial  = getMaterial("toyDS.insideMaterialName");
    
    // Toy model of the DS coils + cryostat. It needs more structure and has
    // much less total material.
    VolumeInfo detSolCoilInfo = nestTubs( "ToyDSCoil",
					  detSolCoilParams,
					  detSolCoilMaterial,
					  0,
					  G4ThreeVector(),
					  shieldFeInsideInfo.logical,
					  0,
					  G4Color::Red()
					  );
    
    // The vaccum inside the DS cryostat and coils; this is longer in z than the coils+cryo.
    VolumeInfo detSolVacInfo = nestTubs( "ToyDSVacuum",
					 detSolVacParams,
					 detSolVacMaterial,
					 0,
					 G4ThreeVector(),
					 shieldFeInsideInfo.logical,
					 0,
					 G4Color::Magenta()
					 );

    // Mock up of the production solenoid and its vacuum.

    double prodSolCoilParams[5] = { 
      _config->getDouble("toyPS.rIn"       ) * mm,
      _config->getDouble("toyPS.rOut"      ) * mm,
      _config->getDouble("toyPS.halfLength") * mm,
      0.,
      2.*M_PI
    };
    double prodSolVacParams[5] = { 
      0.*mm,
      prodSolCoilParams[0],
      _config->getDouble("toyPS.halfLengthVac") * mm,
      0.,
      2.*M_PI
    };
    
    G4Material* prodSolCoilMaterial = getMaterial("toyPS.materialName");
    G4Material* prodSolVacMaterial  = getMaterial("toyPS.insideMaterialName");
    
    // Position of PS inside the air volume of the hall.
    G4ThreeVector prodSolCoilOffset = 
      G4ThreeVector( prodSolXoff-hallPosition[0],
		     yOriginHeight - hallInHLen[1],
		     _config->getDouble("toyPS.z0") * mm + _mu2eOrigin.z()
		     );

    // Toy model of the PS coils + cryostat. It needs more structure and has
    // much less total material.
    VolumeInfo prodSolCoilInfo = nestTubs( "ToyPSCoil",
					   prodSolCoilParams,
					   prodSolCoilMaterial,
					   0,
					   prodSolCoilOffset,
					   hallInfo.logical,
					   0,
					   G4Color::Cyan()
					   );
    
    // The vacuum inside the PS cyrostat; this can be longer than the coils!
    VolumeInfo prodSolVacInfo = nestTubs( "ToyPSVacuum",
					  prodSolVacParams,
					  prodSolVacMaterial,
					  0,
					  prodSolCoilOffset,
					  hallInfo.logical,
					  0,
					  G4Color::Yellow()
					  );

    VolumeInfo trackerInfo;
    //trackerInfo.solid   = new G4Tubs( name, param[0], param[1], param[2], param[3], param[4]  );
    //trackerInfo.logical = new G4LogicalVolume( info.solid, material, name); 


    if( _config->getBool("hasLTracker",false) ){
      trackerInfo = constructLTracker( detSolVacInfo.logical, dsz0 );

    } else{


      // Make a TUBs to represent the tracking volume.
      // Add detail later.
      // A hack here - should get numbers via the geometry system.
      double trackerParams[5] = { 
	0.,
	800.,
	1300.,
	0.,
	2.*M_PI
      };

      G4Material*    trackerMaterial = detSolVacMaterial;
      G4ThreeVector trackerOffset(0.,0.,12000.-dsz0-1800.);
      trackerInfo = nestTubs( "TrackerMother",
			      trackerParams,
			      trackerMaterial,
			      0,
			      trackerOffset,
			      detSolVacInfo.logical,
			      0,
			      G4Color::Yellow(),
			      true
			      );
    }

    // Make a TUBs to represent the target system.
    // Add detail later.
    // A hack here - should get numbers via the geometry system.
    double targetParams[5] = { 
      0.,
      100.,
      400.,
      0.,
      2.*M_PI
    };
    G4Material*   targetMaterial = detSolVacMaterial;
    G4ThreeVector targetOffset(0.,0.,12000.-dsz0-6100.);
    VolumeInfo targetInfo = nestTubs( "TargetMother",
				       targetParams,
				       targetMaterial,
				       0,
				       targetOffset,
				       detSolVacInfo.logical,
				       0,
				       G4Color::Yellow(),
				       true
				       );

    // Only after all volumes have been defined should we set the magnetic fields.
    // (Or else we need to set it on a volume by volume basis; this way we can
    //  set it on one volume and all of its subvolumes.)

    // Make the magnetic field valid inside the detSol vacuum.
    G4double bz = _config->getDouble("toyDS.bz") * tesla;
    _detSolBField = auto_ptr<G4UniformMagField>(new G4UniformMagField(G4ThreeVector(0.,0.,bz)));
    
    // Create a field manager to manange this field.  Use exact helix stepping.
    G4double stepMinimum(1.0e-2*mm);
    _usualRHS    = auto_ptr<G4Mag_UsualEqRhs>   (new G4Mag_UsualEqRhs( _detSolBField.get()) );
    _exactHelix  = auto_ptr<G4ExactHelixStepper>(new G4ExactHelixStepper(_usualRHS.get()));
    _chordFinder = auto_ptr<G4ChordFinder>      (new G4ChordFinder( _detSolBField.get(), stepMinimum, _exactHelix.get() ));
    _fieldMgr    = auto_ptr<G4FieldManager>     (new G4FieldManager( _detSolBField.get(), _chordFinder.get(), true));

    // Attach the field manager to the detSol volume.
    detSolVacInfo.logical->SetFieldManager( _fieldMgr.get(), true);
    
    // Set step limit.  
    // See also PhysicsList.cc to add a steplimiter to the list of processes.
    // Do this so that we can see the helical trajectory in the DS and volumes inside of it.
    G4double maxStep = 20.;
    _stepLimit = auto_ptr<G4UserLimits>( new G4UserLimits(maxStep));
    detSolVacInfo.logical->SetUserLimits(_stepLimit.get());
    trackerInfo.logical->SetUserLimits(_stepLimit.get());
    targetInfo.logical->SetUserLimits(_stepLimit.get());

  }

  // Get the requested material, throwing if necessary.
  G4Material* Mu2eWorld::getMaterial( string const& key ){
    string materialName = _config->getString(key);
    return findMaterialOrThrow(materialName);
  }

  // Convert to base units for all of the items in the vector.
  void Mu2eWorld::setUnits( vector<double>& V, G4double unit ){
    for ( vector<double>::iterator b=V.begin(), e=V.end();
	  b!=e; ++b){
      *b *= unit;
    }
  }

  // Helper function used to construct the LTracker.
  void printStrawInfo( Straw const& s){
    cout << "Straw: " 
	 << s.getMidPoint()
	 << endl;
  }

  VolumeInfo Mu2eWorld::constructLTracker( G4LogicalVolume* mother, double zOff ){

    // Master geometry for the LTracker.
    GeomHandle<LTracker> ltracker;

    double rOut  = mm * ltracker->rOut();
    double zHalf = mm * ltracker->zHalfLength();
    double z0    = mm * ltracker->z0();

    VolumeInfo trackerInfo;

    // Make the mother volume for the LTracker.
    string trackerName("LTrackerMother");
    G4Material* fillMaterial = findMaterialOrThrow(ltracker->fillMaterial());
    G4ThreeVector trackerOffset(0.,0.,z0-zOff);

    trackerInfo.solid  = new G4Tubs( trackerName,
				     0., rOut, zHalf, 0., 2.*M_PI );
    
    trackerInfo.logical = new G4LogicalVolume( trackerInfo.solid, fillMaterial, trackerName); 
    
    trackerInfo.physical =  new G4PVPlacement( 0, 
					       trackerOffset, 
					       trackerInfo.logical, 
					       trackerName, 
					       mother, 
					       0, 
					       0);

    // Visualization attributes of the the mother volume.
    {
      G4VisAttributes* visAtt = new G4VisAttributes(true, G4Colour::Green() );
      visAtt->SetForceSolid(true);
      visAtt->SetForceAuxEdgeVisible (false);
      visAtt->SetVisibility(true);
      trackerInfo.logical->SetVisAttributes(visAtt);
    }

    // For now cheat and assume that all straws are the same and are just made of gas with
    // no walls or wires.
    Straw const& straw = ltracker->getStraw( StrawId( LTracker::wedge, 0, 0, 0) );
    StrawDetail const& detail = straw.getDetail();

    string strawName("Straw");
    VolumeInfo strawInfo;
      
    G4Material* strawMaterial = findMaterialOrThrow( detail.materialName(1) );
    strawInfo.solid  = new G4Tubs(strawName
				  ,0.
				  ,detail.outerRadius() * mm
				  ,detail.halfLength()  * mm
				  ,0.
				  ,CLHEP::twopi*radian
				  );
    
    strawInfo.logical = new G4LogicalVolume( strawInfo.solid
					     , strawMaterial
					     , strawName
					     );

    // Define the straws to be sensitive detectors.
    // Does this leak the SDman?
    G4SDManager* SDman   = G4SDManager::GetSDMpointer();
    G4String strawSDname = "StrawGasVolume";
    StrawSD* strawSD     = new StrawSD( strawSDname );
    SDman->AddNewDetector( strawSD );
    strawInfo.logical->SetSensitiveDetector( strawSD );


    // Does this leak strawVisAtt ??
    {
      G4VisAttributes* strawVisAtt = new G4VisAttributes(true, G4Colour::Green() );
      strawVisAtt->SetForceSolid(false);
      strawVisAtt->SetForceAuxEdgeVisible (false);
      strawVisAtt->SetVisibility(false);
      strawInfo.logical->SetVisAttributes(strawVisAtt);
    }

    // Cheat again and do not bother to segment the tracker volume.
    // Just place the straws in the final position.
    // Need to properly segment the volume on a per sector basis in
    // order to improve G4 statup speed.
    StrawPlacer placer( "StrawPhys", strawInfo.logical, trackerInfo.logical );
    ltracker->forAllStraws( placer);

    return trackerInfo;

  }


  
} // end namespace mu2e
