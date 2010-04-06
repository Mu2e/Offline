//
// Construct the Mu2e G4 world and serve information about that world.
//
// $Id: Mu2eWorld.cc,v 1.14 2010/04/06 19:27:19 rhbob Exp $
// $Author: rhbob $ 
// $Date: 2010/04/06 19:27:19 $
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

#include <boost/regex.hpp>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/StrawPlacer.hh"
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/ITGasLayerSD_ExtWireData.hh"
#include "Mu2eG4/inc/ITGasLayerSD_v2.hh"
#include "Mu2eG4/inc/ITGasLayerSD_v3.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/ITrackerBuilder.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "TargetGeom/inc/Target.hh"


// G4 includes
#include "G4AssemblyVolume.hh"
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
#include "G4Hype.hh"
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
#include "G4ClassicalRK4.hh"
#include "G4CashKarpRKF45.hh"
#include "G4ImplicitEuler.hh"
#include "G4ExplicitEuler.hh"

//#include "G4GDMLParser.hh"

//
//Julie's DSField for reading in the file for full DS field
#include "Mu2eG4/inc/DSField.hh"


using namespace std;

namespace mu2e {

  // Value used to request that all Mu2e specific materials be made.
  static const std::string DoAllValue = "DoAll";

  Mu2eWorld::Mu2eWorld()
    :  _cosmicReferencePoint(),
       _mu2eOrigin(),
       _info(),
       _detSolUpstreamBField(),
       _usualUpstreamRHS(),
       _exactUpstreamHelix(),
       _chordUpstreamFinder(),
       _fieldUpstreamMgr(),
       _stepUpstreamLimit(),
       _detSolDownstreamBField(),
       _usualDownstreamRHS(),
       _exactDownstreamHelix(),
       _chordDownstreamFinder(),
       _fieldDownstreamMgr(),
       _stepDownstreamLimit(),
       _lTrackerWedgeAssembly(),
       _lTrackerVaneAssembly(),
       _lTrackerAssemblyVols(){
  }
  
  Mu2eWorld::~Mu2eWorld(){

    // Do not destruct the solids, logical volumes or physical volumes.
    // G4 looks after that itself.

  }

  void printPhys() {
    G4PhysicalVolumeStore* pstore = G4PhysicalVolumeStore::GetInstance();
    int n(0);
    for ( std::vector<G4VPhysicalVolume*>::const_iterator i=pstore->begin(); i!=pstore->end(); i++){
      cout << "Physical Volume: "
	   << setw(5) << n++
	   << (*i)->GetName()
	   << endl;
      if ( n > 25 ) break;
    }

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
    edm::LogInfo log("GEOM");
    log << "Mu2e Origin: " << _mu2eOrigin << "\n";

    // Origin used to construct the MECO detector.
    _mu2eDetectorOrigin = _mu2eOrigin + G4ThreeVector( -3904., 0., 12000.);
    log << "Mu2e Detector Origin: " << _mu2eDetectorOrigin << "\n";

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
    log << "Cosmic Ref: " << _cosmicReferencePoint << "\n";

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
    
    // this was supposed to be in/out,  but it looks like toyDS.rOut and rIn are reversed, 
    //so I will put them in the wrong order.  kutschke agrees was a bug.
     double detSolCoilParams[5] = { 
      _config->getDouble("toyDS.rIn"       ) * mm,
      _config->getDouble("toyDS.rOut"      ) * mm,
      _config->getDouble("toyDS.halfLength") * mm,
      0.,
      2.*M_PI
    };
 
   //will divide the detector solenoid vacuum into two parts. the purpose
    //is to have a slowly falling field (the real one) over the stopping target
    //region, and then a pure solenoidal field in the region of the tracker.
    //
    //this will greatly simplify debugging of the kalman filter and tracking algorithms.
    //I considered fudging the field Julie Managan put in, which is the upstream half; Kutschke
    //says this will cause problems because no matter what I do it will violate Poisson's equation
    //at some point.  By dividing it into two parts, I can attach different fieldMgrs to the two sections
    //and this does less emotional violence to G4 as long as I choose a nice transition point.
    // rhb 1/22/10
    //
    //this transitionZ is in theory the only hardwired fudge to make this work; then calculate everything
    //given where we choose to move from one fieldMgr to the other

    //
    //ok, now compute the new halfLengths and centers; note Rob's code forces center at zero locally but give it a name for readability
    double centerOfDS = 0.;
    //now, you can read off a transition Z based on the field map and translate it to the local system
    double globalTransitionZ = dsz0; //so let's start at local 0 for debugging
    double transitionZ = (globalTransitionZ - dsz0) - 1500.;//transition arbitrary while I debug but must respect target/tracker locations
    
    double halfLengthOfUpstreamDSVac     = (_config->getDouble("toyDS.halfLengthVac")*mm + (transitionZ-centerOfDS))/2.;
    double halfLengthOfDownstreamDSVac   = _config->getDouble("toyDS.halfLengthVac")*mm - halfLengthOfUpstreamDSVac;
    //and of course the center of the upstream and downstream sections just moved
    double centerOfUpstreamDSVac   = transitionZ - halfLengthOfUpstreamDSVac;
    double centerOfDownstreamDSVac = transitionZ + halfLengthOfDownstreamDSVac;

    //print out all this nonsense
    /*
    cout << "we started with a center at: " << globalTransitionZ - dsz0
	 << " and a halflength of " << _config->getDouble("toyDS.halfLengthVac") 
	 << " a transition z at " << transitionZ 
	 << " and ended up with " << endl;

    cout << " halfLengthOfUpstreamDSVac = " << halfLengthOfUpstreamDSVac
	 << " centerOfUpstreamDSVac = " << centerOfUpstreamDSVac
	 << endl;
    cout << " halfLengthOfDownstreamDSVac = " << halfLengthOfDownstreamDSVac
	 << " centerOfDownstreamDSVac = " << centerOfDownstreamDSVac
	 << endl;
    */
    //checked this out for a special case, looked fine
    //     assert(2==1);

    double detSolVacParams[5] = { 
      0. * mm,
      detSolCoilParams[0],
      _config->getDouble("toyDS.halfLengthVac") * mm,
      0.,
      2.*M_PI
      };

    //cout << "toyDS.halfLengthVac = " << _config->getDouble("toyDS.halfLengthVac") << endl;
    //    assert (2==1);
    double detSolUpstreamVacParams[5]   = { 
       0.
      ,detSolCoilParams[0]
      ,halfLengthOfUpstreamDSVac*mm
      ,0.
      ,2.*M_PI
    };
    double detSolDownstreamVacParams[5]   = { 
       0.
      ,detSolCoilParams[0]
      ,halfLengthOfDownstreamDSVac*mm
      ,0.
      ,2.*M_PI
    };
    G4Material* detSolCoilMaterial = getMaterial("toyDS.materialName");
    //no longer used    G4Material* detSolVacMaterial  = getMaterial("toyDS.insideMaterialName");
    G4Material* detSolUpstreamVacMaterial    = getMaterial("toyDS.insideMaterialName");
    G4Material* detSolDownstreamVacMaterial  = getMaterial("toyDS.insideMaterialName");
    
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
    /* again not used with split 
    // The vacuum inside the DS cryostat and coils; this is longer in z than the coils+cryo.
    VolumeInfo detSolVacInfo = nestTubs( "ToyDSVacuum",
					 detSolVacParams,
					 detSolVacMaterial,
					 0,
					 G4ThreeVector(),
					 shieldFeInsideInfo.logical,
					 0,
					 G4Color::Magenta()
					 );
    */


    G4ThreeVector detSolDownstreamOffset  = G4ThreeVector(0.,0.,centerOfDownstreamDSVac);
    VolumeInfo detSolDownstreamVacInfo = nestTubs( "ToyDSDownstreamVacuum",
					 detSolDownstreamVacParams,
					 detSolDownstreamVacMaterial,
					 0,
					 detSolDownstreamOffset,
					 shieldFeInsideInfo.logical,
					 0,
					 G4Color::Magenta()
					 );

    G4ThreeVector detSolUpstreamOffset  = G4ThreeVector(0.,0.,centerOfUpstreamDSVac);
    VolumeInfo detSolUpstreamVacInfo   = nestTubs( "ToyDSUpstreamVacuum",
					 detSolUpstreamVacParams,
					 detSolUpstreamVacMaterial,
					 0,
					 detSolUpstreamOffset,
					 shieldFeInsideInfo.logical,
					 0,
					 G4Colour::White() //color change between two halves
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


   // Proton Target in PS 

    // Proton Target parameters 
    // Proton Target position
    G4ThreeVector ProtonTargetPosition = G4ThreeVector( 
				_config->getDouble("targetPS_positionX")*mm,
				_config->getDouble("targetPS_positionY")*mm,
				_config->getDouble("targetPS_positionZ")*mm
     				);
    
    // Rotation of Proton Target				
    double targetPS_rotX = _config->getDouble("targetPS_rotX" );
    double targetPS_rotY = _config->getDouble("targetPS_rotY" );
    
    // Proton Target Material
    G4Material* targetPS_materialName = getMaterial("targetPS_materialName");   
    
    //Proton Target geometry parameters
    double targetPS_Pam[5] = { 
      0.,
      _config->getDouble("targetPS_rOut"      ) * mm,
      _config->getDouble("targetPS_halfLength") * mm,
      0.,
      2.*M_PI
    };

   // Rotation of Proton Target
   G4RotationMatrix* PS_target_rot = new G4RotationMatrix();
    PS_target_rot->rotateX( targetPS_rotX*degree);
    PS_target_rot->rotateY( targetPS_rotY*degree);
   
    //
    // Creating Proton Target object in prodSolVacInfo logical object (v. khalatian)
   VolumeInfo ProtonTargetInfo = nestTubs( "ProtonTarget",
					  targetPS_Pam,
					  targetPS_materialName,
					  PS_target_rot,
					  ProtonTargetPosition,
					  prodSolVacInfo.logical,
					  0,
					  G4Color::White()
					  );
   
    // Primary Proton Gun Origin 
    _primaryProtonGunOrigin = dirtOffset + wallOffset + hallOffset + prodSolCoilOffset + ProtonTargetPosition;
    
    //Primary Proton Gun Rotation 
    // For rotating Primary Proton Gun I take angles from Proton Target 
    _primaryProtonGunRotation.rotateX( targetPS_rotX*degree);
    _primaryProtonGunRotation.rotateY( targetPS_rotY*degree);
    //
    //these are "active rotations; we want passive, in G4 style
    _primaryProtonGunRotation = _primaryProtonGunRotation.inverse();


    VolumeInfo trackerInfo;
    //trackerInfo.solid   = new G4Tubs( name, param[0], param[1], param[2], param[3], param[4]  );
    //trackerInfo.logical = new G4LogicalVolume( info.solid, material, name); 

    if( _config->getBool("hasLTracker",false) ){
      int ver = _config->getInt("LTrackerVersion",3);
      log << "LTracker version: " << ver << "\n";

// kutschke says use only v3, significant startup and performance penalty 
// for v2. default was set to 2
// beign careful to associate with downstream detsol
      log << "LTracker version: " << ver << "\n";
      if ( ver == 0 ){
	trackerInfo = constructLTracker( detSolDownstreamVacInfo.logical, dsz0 + centerOfDownstreamDSVac );
      }
      else if ( ver == 1 ) {
	trackerInfo = constructLTrackerv2( detSolDownstreamVacInfo.logical, dsz0 + centerOfDownstreamDSVac );
      } else {
	trackerInfo = constructLTrackerv3( detSolDownstreamVacInfo.logical, dsz0 + centerOfDownstreamDSVac );
      }
    } else if ( _config->getBool("hasITracker",false) ) {
    	trackerInfo = ITrackerBuilder::constructTracker( detSolDownstreamVacInfo.logical, dsz0 + centerOfDownstreamDSVac );
    } else {

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

    	//
    	//tracker is associated with downstream constant section -- that's the point
    	G4Material*    trackerMaterial = detSolDownstreamVacMaterial;

    	//
    	// this doesn't need to change -- the tracker doesn't move.
    	// however, it really should be computed

    	double trackerCenterInZ = 12000. - dsz0 - 1800.;
    	G4ThreeVector trackerOffset(0.,0.,trackerCenterInZ);
    	/*
           cout << "tracker center in Z " << trackerCenterInZ << "\n" 
	   <<trackerParams[0] << "\n " 
	   <<trackerParams[1] << "\n " 
	   <<trackerParams[2] << "\n" 
	   <<trackerParams[3] << "\n" 
	   <<trackerParams[4] << "\n" 
	   <<endl;
    	 */
    	if ( (trackerCenterInZ - trackerParams[2]) < transitionZ)
    	{ cout << "from Mu2eWorld:  transition Z in the middle of the tracker, fool..." << endl;
    	assert(2==1);
    	}

    	trackerInfo = nestTubs( "TrackerMother",
    			trackerParams,
    			trackerMaterial,
    			0,
    			trackerOffset,
    			detSolDownstreamVacInfo.logical,
    			0,
    			G4Color::Yellow(),
    			true
    	);


    }

    // Do the Target

    VolumeInfo targetInfo;

    if( _config->getBool("hasTarget",false) ){

	targetInfo = constructTarget( detSolUpstreamVacInfo.logical, 
                        dsz0 + centerOfUpstreamDSVac );

    } else {
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

    // 12000 - dsz0 - 6100 is from kutschke when the mother of the target was detSolVacInfo;
    double targetCenterInZ = 12000. - dsz0 - 6100. - centerOfUpstreamDSVac;
    G4ThreeVector targetOffset(0.,0.,targetCenterInZ);
    //make sure transition  between upstream/downstream isn't in the target since
    //the target is only in the upstream part
      cout << "target center in Z= " << targetCenterInZ<< "\n" 
	   <<targetParams[0] << "\n" 
	   <<targetParams[1] << "\n" 
	   <<targetParams[2] << "\n"
	   <<targetParams[3] << "\n" 
	   <<targetParams[4] << "\n" 
	   <<endl;

      if ( (targetCenterInZ + targetParams[2])+centerOfUpstreamDSVac > transitionZ)
	{ cout << "from Mu2eWorld:  transition Z in the middle of the target, idiot (Bernstein!)..." << endl;
	  assert(false);
	}

      G4Material*   targetMaterial = detSolUpstreamVacMaterial;
      targetInfo = nestTubs( "TargetMother",
				       targetParams,
				       targetMaterial,
				       0,
				       targetOffset,
				       detSolUpstreamVacInfo.logical,
				       0,
				       G4Color::Yellow(),
				       true
				       );
    } //hasTarget

    // Target done.

    // Only after all volumes have been defined should we set the magnetic fields.
    // Make the magnetic field valid inside the detSol vacuum; one upstream, one downstream

    const char* fieldmap = "/home2/misc1/jmanagan/myMu2e/GMC/fieldmaps/dsmap_unfmt_rad100.dat";//note disgusting hardwired absolute path
    int const nx(50); int const ny(25); int const nz(438);

    G4double stepUpstreamMinimum(1.0e-2*mm);
    G4double stepDownstreamMinimum(1.0e-2*mm);

    //
    //get the field form; default if unspecified is constant
    int detSolFieldForm = _config->getInt("detSolFieldForm",detSolUpConstantDownConstant); 

      cout << "detSolFieldForm  from Mu2eWorld.cc = " << detSolFieldForm << endl;
    //assert (2==1);
    // first check we have a legal configuration
    if (detSolFieldForm != detSolFullField && detSolFieldForm != detSolUpVaryingDownConstant && detSolFieldForm != detSolUpConstantDownConstant)
      {
	G4cout << " no legal field specification; detSolFieldForm = " << detSolFieldForm << G4endl;
      throw cms::Exception("GEOM")
	<< "illegal field config as specified in geom.txt \n";
      }

    if (detSolFieldForm == detSolFullField)
      {
	//
	//upstream varying section

	_detSolUpstreamVaryingBField = auto_ptr<DSField>(new DSField(fieldmap,_mu2eOrigin,nx,ny,nz));
	_usualUpstreamRHS    = auto_ptr<G4Mag_UsualEqRhs>   (new G4Mag_UsualEqRhs( _detSolUpstreamVaryingBField.get() ) );
	_rungeEEUpstreamHelix  = auto_ptr<G4ExplicitEuler>(new G4ExplicitEuler(_usualUpstreamRHS.get()));
        _chordUpstreamFinder = auto_ptr<G4ChordFinder>      (new G4ChordFinder( _detSolUpstreamVaryingBField.get(), stepUpstreamMinimum
									    , _rungeUpstreamHelix.get() ));
        _fieldUpstreamMgr    = auto_ptr<G4FieldManager>     (new G4FieldManager( _detSolUpstreamVaryingBField.get(), 
										 _chordUpstreamFinder.get(), true));
	//
	//downstream varying section
	_detSolDownstreamVaryingBField = auto_ptr<DSField>(new DSField(fieldmap,_mu2eOrigin,nx,ny,nz));
	_usualDownstreamRHS    = auto_ptr<G4Mag_UsualEqRhs>   (new G4Mag_UsualEqRhs( _detSolDownstreamVaryingBField.get() ) );
	_rungeEEDownstreamHelix  = auto_ptr<G4ExplicitEuler>(new G4ExplicitEuler(_usualDownstreamRHS.get()));
        _chordDownstreamFinder = auto_ptr<G4ChordFinder>      (new G4ChordFinder( _detSolDownstreamVaryingBField.get(), stepDownstreamMinimum,
										  _rungeDownstreamHelix.get() ));
        _fieldDownstreamMgr    = auto_ptr<G4FieldManager>     (new G4FieldManager( _detSolDownstreamVaryingBField.get(), 
										   _chordDownstreamFinder.get(), true));

      }
    if (detSolFieldForm == detSolUpVaryingDownConstant)
      {
	cout << "in hybrid " << endl;
	//
	//upstream varying section
	_detSolUpstreamVaryingBField = auto_ptr<DSField>(new DSField(fieldmap,_mu2eOrigin,nx,ny,nz));
	_usualUpstreamRHS    = auto_ptr<G4Mag_UsualEqRhs>   (new G4Mag_UsualEqRhs( _detSolUpstreamVaryingBField.get()));
	
	_rungeEEUpstreamHelix  = auto_ptr<G4ExplicitEuler> (new G4ExplicitEuler(_usualUpstreamRHS.get()));
        _chordUpstreamFinder = auto_ptr<G4ChordFinder>      (new G4ChordFinder( _detSolUpstreamVaryingBField.get(), stepUpstreamMinimum
									    , _rungeEEUpstreamHelix.get() ));
        _fieldUpstreamMgr    = auto_ptr<G4FieldManager>     (new G4FieldManager( _detSolUpstreamVaryingBField.get(), 
								 _chordUpstreamFinder.get(), true));

	//downstream constant section
	G4double bzDown = _config->getDouble("toyDS.bz") * tesla;
	_detSolDownstreamConstantBField = auto_ptr<G4UniformMagField>(new G4UniformMagField(G4ThreeVector(0.,0.,bzDown)));
	_usualDownstreamRHS    = auto_ptr<G4Mag_UsualEqRhs>   (new G4Mag_UsualEqRhs( _detSolDownstreamConstantBField.get()) );
	_exactDownstreamHelix  = auto_ptr<G4ExactHelixStepper>(new G4ExactHelixStepper(_usualDownstreamRHS.get()));
        _chordDownstreamFinder = auto_ptr<G4ChordFinder>      (new G4ChordFinder( _detSolDownstreamConstantBField.get(), stepDownstreamMinimum
									      , _exactDownstreamHelix.get() ));
        _fieldDownstreamMgr    = auto_ptr<G4FieldManager>     (new G4FieldManager( _detSolDownstreamConstantBField.get(), 
										   _chordDownstreamFinder.get(), true));
      }

    if (detSolFieldForm == detSolUpConstantDownConstant) 
      {
	cout << "in constant field" << endl;
	//
	// constant field, but split into two parts; upstream first
	G4double bzUp = _config->getDouble("toyDS.bz") * tesla;
	_detSolUpstreamConstantBField = auto_ptr<G4UniformMagField>(new G4UniformMagField(G4ThreeVector(0.,0.,bzUp)));
	_usualUpstreamRHS    = auto_ptr<G4Mag_UsualEqRhs>   (new G4Mag_UsualEqRhs( _detSolUpstreamConstantBField.get() ) );
	_exactUpstreamHelix  = auto_ptr<G4ExactHelixStepper>(new G4ExactHelixStepper(_usualUpstreamRHS.get()));
        _chordUpstreamFinder = auto_ptr<G4ChordFinder>      (new G4ChordFinder( _detSolUpstreamConstantBField.get(), stepUpstreamMinimum
										, _exactUpstreamHelix.get() ));
        _fieldUpstreamMgr    = auto_ptr<G4FieldManager>     (new G4FieldManager( _detSolUpstreamConstantBField.get(), 
										 _chordUpstreamFinder.get(), true));
	//downstream
	G4double bzDown = _config->getDouble("toyDS.bz") * tesla;
	_detSolDownstreamConstantBField = auto_ptr<G4UniformMagField>(new G4UniformMagField(G4ThreeVector(0.,0.,bzDown)));
	_usualDownstreamRHS    = auto_ptr<G4Mag_UsualEqRhs>   (new G4Mag_UsualEqRhs( _detSolDownstreamConstantBField.get() ) );
	_usualDownstreamRHS    = auto_ptr<G4Mag_UsualEqRhs>   (new G4Mag_UsualEqRhs( _detSolDownstreamConstantBField.get() ) );
	_exactDownstreamHelix  = auto_ptr<G4ExactHelixStepper>(new G4ExactHelixStepper(_usualDownstreamRHS.get()));
        _chordDownstreamFinder = auto_ptr<G4ChordFinder>      (new G4ChordFinder( _detSolDownstreamConstantBField.get(), stepDownstreamMinimum
										, _exactDownstreamHelix.get() ));
        _fieldDownstreamMgr    = auto_ptr<G4FieldManager>     (new G4FieldManager( _detSolDownstreamConstantBField.get(), 
										   _chordDownstreamFinder.get(), true));
      }



    // Now that we've chosen, attach the field manager to the detSol volume; full field upstream
    detSolUpstreamVacInfo.logical->SetFieldManager( _fieldUpstreamMgr.get(), true);
    detSolDownstreamVacInfo.logical->SetFieldManager( _fieldDownstreamMgr.get(), true);



    //
    //set integration step values
    G4double singleValue = 0.5e-02*mm;
    G4double newUpstreamDeltaI = singleValue;
    G4double newDownstreamDeltaI = singleValue;
    G4double deltaOneStep = singleValue;
    G4double deltaChord = singleValue;
    G4double maxStep = 1.e-00*mm;

    _fieldUpstreamMgr->SetDeltaOneStep(deltaOneStep);
    _fieldDownstreamMgr->SetDeltaOneStep(deltaOneStep);
    

    _fieldUpstreamMgr->SetDeltaIntersection(newUpstreamDeltaI);    
    _fieldDownstreamMgr->SetDeltaIntersection(newDownstreamDeltaI);    

    _chordUpstreamFinder->SetDeltaChord(deltaChord);
    _chordDownstreamFinder->SetDeltaChord(deltaChord);


    
    // Set step limit.  
    // See also PhysicsList.cc to add a steplimiter to the list of processes.
    // Do this so that we can see the helical trajectory in the DS and volumes inside of it.
    _stepLimit = auto_ptr<G4UserLimits>( new G4UserLimits(maxStep));
    detSolUpstreamVacInfo.logical->SetUserLimits(_stepLimit.get());
    detSolDownstreamVacInfo.logical->SetUserLimits(_stepLimit.get());

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

  // Option 1:
  // Build LTracker with no substructure.  
  // Place each straw in the tracker mother volume.
  //
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

  // Option 2:
  // Build LTracker from assembly volumes.
  // 1) Make one assembly volume for the octant sides and one for the vanes.
  // 2) Use imprint to make 8 copies of each and put each in the right location.
  //
  VolumeInfo Mu2eWorld::constructLTrackerv2( G4LogicalVolume* mother, double zOff ){
    //    VolumeInfo trackerInfo;

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


    /*
    cout << "Tracker Offset: z0, zOff, z0-zOff: " 
	 << z0 << " "
	 << zOff << " "
	 << z0-zOff << " "
	 << endl;
    */

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

    Straw const& straw        = ltracker->getStraw( StrawId( LTracker::wedge, 0, 0, 0) );
    StrawDetail const& detail = straw.getDetail();

    // Build logical volume for a straw.
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

    if ( ltracker->getDevices().size() != ndevices ){
      throw cms::Exception("GEOM")
        << "Unexpected number of devices in the LTracker: "
        << ltracker->getDevices().size()
        << "\n";
    }

    // Build an assembly volume for each of the octants and vanes.
    for ( std::size_t idev=0; idev<ltracker->getDevices().size(); ++idev ){

      // Assume all sectors are the same as sector 0.
      Sector const& sec = ltracker->getSector(SectorId(idev,0));
	
      _lTrackerAssemblyVols[idev] = auto_ptr<G4AssemblyVolume> (new G4AssemblyVolume());

      for ( std::size_t ilay =0; ilay<sec.getLayers().size(); ++ilay){
	Layer const& lay = sec.getLayer(ilay);
	Hep3Vector const& origin = sec.getBasePosition().at(ilay);
	Hep3Vector const& delta  = sec.getBaseDelta();

	StrawId id(idev,0,ilay,0);
	
	for ( std::size_t istr = 0; istr<lay.nStraws(); ++istr ){
	  Hep3Vector position = origin + istr*delta;
	  _lTrackerAssemblyVols[idev]->AddPlacedVolume( strawInfo.logical, position, 0); 
	}
      }
      
    }

    // Imprint the assembly volumes.
    for ( std::size_t idev = 0; idev<ltracker->getDevices().size(); ++idev){
      Device const& device = ltracker->getDevice(idev);

      for ( std::size_t isec =0; isec<device.getSectors().size(); ++isec){
	Sector const& sector = device.getSector(isec);

	HepRotationX RX(sector.boxRxAngle());
	HepRotationY RY(sector.boxRyAngle());
	HepRotationZ RZ(sector.boxRzAngle());
     
	// Need to understand if this causes memory leak.
	G4RotationMatrix* rot = new G4RotationMatrix( RZ*RX*RY);

	// MakeImprint requires non-const argument.
	G4ThreeVector offset = sector.boxOffset();


	// Copy numbers start at base+1
	StrawId id(idev,isec,0,0);
	int baseCopyNumber = ltracker->getStraw(id).Index().asInt()-1;

	_lTrackerAssemblyVols[idev]->MakeImprint( trackerInfo.logical, offset, rot, baseCopyNumber); 

      }
    }

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
    
    return trackerInfo;

  }

  // Version 3 of LTracker.
  // Make boxes to bound each sector.
  // Place straws within each box.
  VolumeInfo Mu2eWorld::constructLTrackerv3( G4LogicalVolume* mother, double zOff ){

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

    /*
    cout << "Tracker Offset: z0, zOff, z0-zOff: " 
	 << z0 << " "
	 << zOff << " "
	 << z0-zOff << " "
	 << endl;
    */

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

    Straw const& straw        = ltracker->getStraw( StrawId( LTracker::wedge, 0, 0, 0) );
    StrawDetail const& detail = straw.getDetail();

    // Build logical volume for a straw.
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

    vector<VolumeInfo> vinfo;

    for ( std::size_t idev = 0; idev<ltracker->getDevices().size(); ++idev){
      Device const& device = ltracker->getDevice(idev);

      for ( std::size_t isec =0; isec<device.getSectors().size(); ++isec){
	Sector const& sector = device.getSector(isec);

	// Name of this sector as string.
	string name = sector.name("LTrackerSector_");

	// Construct the rotation.  
	// This rotation is the inverse of the one in v2.
	// Note the sign and the reversed order : active/passive  confusion.
	// Need to understand if this causes memory leak.
	HepRotationX RX(-sector.boxRxAngle());
	HepRotationY RY(-sector.boxRyAngle());
	HepRotationZ RZ(-sector.boxRzAngle());
	G4RotationMatrix* rot = new G4RotationMatrix( RY*RX*RZ);

	// Make a physical volume for this sector.  Same material as the 
	// main LTracker volume ( some sort of vacuum ).
	VolumeInfo tmp = nestBox( name,
				  sector.boxHalfLengths(),
				  fillMaterial,
				  rot,
				  sector.boxOffset(),
				  trackerInfo.logical,
				  0);
	vinfo.push_back(tmp);
	VolumeInfo const& sectorBoxInfo = vinfo.back();

	Hep3Vector const& delta  = sector.getBaseDelta();

	for ( std::size_t ilay =0; ilay<sector.getLayers().size(); ++ilay){
	  Layer const& layer = sector.getLayer(ilay);

	  Hep3Vector const& origin = sector.getBasePosition().at(ilay);

	  for ( std::size_t istr =0; istr<layer.nStraws(); ++istr){
	    Straw const& straw = layer.getStraw(istr);

	    // Position within the sector box.
	    Hep3Vector position = origin + istr*delta;

	    // Name of this physical volume.
	    string sname = straw.name( "LTrackerStraw_");

	    G4VPhysicalVolume* phys = new G4PVPlacement( 0, 
							 position,
							 strawInfo.logical,
							 sname, 
							 sectorBoxInfo.logical, 
							 0, 
							 straw.Index().asInt()
							 );
	    
	  } // loop over straws
	}   // loop over layers
	
      } // loop over sectors
    }   // loop over devices

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
    
    return trackerInfo;

  }

  VolumeInfo Mu2eWorld::constructTarget( G4LogicalVolume* mother, double zOff ){

std::cout<<"In constructTarget"<<std::endl;
    // Master geometry for the Target assembly
    GeomHandle<Target> target;

    double rOut  = mm * target->cylinderRadius();
    double zHalf = mm * target->cylinderLength()/2.;

    // center in detector coords, assumed to be on axis
    double z0    = mm * target->cylinderCenter();

    VolumeInfo targetInfo;

    // Make the mother volume for the Target
    string targetName("TargetMother");
std::cout<<"Looking for material "<<target->fillMaterial()<<std::endl;
    G4Material* fillMaterial = findMaterialOrThrow(target->fillMaterial());
std::cout<<"Done Looking for material "<<target->fillMaterial()<<std::endl;
    G4ThreeVector targetOffset(0.,0.,(12000+z0-zOff));

    cout << "Target Offset: z0, zOff, z0-zOff: " 
	 << z0 << " "
	 << zOff << " "
	 << z0-zOff << " "
	 << endl;
   


    targetInfo.solid  = new G4Tubs( targetName,
				     0., rOut, zHalf, 0., 2.*M_PI );
    
    targetInfo.logical = new G4LogicalVolume( targetInfo.solid, fillMaterial, targetName); 
    
std::cout<<"targetOffset="<<targetOffset<<std::endl;
std::cout<<"mother has "<<mother->GetNoDaughters()<<" daughters"<<std::endl;
std::cout<<" they are:"<<std::endl;
for (int id=0; id<mother->GetNoDaughters(); id++) cout<<id<<"="<<
      mother->GetDaughter(id)->GetName()<<" at "<<mother->GetDaughter(id)->GetTranslation()<<std::endl;
    targetInfo.physical =  new G4PVPlacement( 0, 
					       targetOffset, 
					       targetInfo.logical, 
					       targetName, 
					       mother, 
					       0, 
					       0);

    // Visualization attributes of the the mother volume.
    {
      // i.e., none...
      G4VisAttributes* visAtt = new G4VisAttributes(false);
      targetInfo.logical->SetVisAttributes(visAtt);
    }

    // now create the individual targets

    for (unsigned int itf=0; itf<target->nFoils(); itf++)
    {

        TargetFoil foil=target->foil(itf);
        VolumeInfo foilInfo;
        G4Material* foilMaterial = findMaterialOrThrow( foil.material() );
        string foilName("Foil");

        foilInfo.solid = new G4Tubs(foilName
                                   ,foil.rIn()
                                   ,foil.rOut()
                                   ,foil.halfThickness()
                                   ,0.
                                   ,CLHEP::twopi*radian
                                   );

        foilInfo.logical = new G4LogicalVolume( foilInfo.solid
                                              , foilMaterial
                                              , foilName
                                              );

        // rotation matrix... 
        G4RotationMatrix* rot = 0; //... will have to wait

        G4ThreeVector foilOffset(foil.center()-G4ThreeVector(0.,0.,z0));
 cout<<"foil "<<itf<<" center="<<foil.center()<<", offset="<<foilOffset<<endl;

        G4VPhysicalVolume* phys = new G4PVPlacement( rot
                                                   , foilOffset
                                                   , foilInfo.logical
                                                   , "TargetFoil_"
                                                   , targetInfo.logical
                                                   , 0
                                                   , itf
                                                   );

        G4VisAttributes* visAtt = new G4VisAttributes(true, G4Colour::Magenta() );
        visAtt->SetForceSolid(true);
        visAtt->SetForceAuxEdgeVisible (false);
        visAtt->SetVisibility(true);
        foilInfo.logical->SetVisAttributes(visAtt);
       



    }// target foils

    return targetInfo;

  }
} // end namespace mu2e
