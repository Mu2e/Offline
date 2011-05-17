//
// Free function to create  Production Solenoid and Production Target.
//
// $Id: constructPS.cc,v 1.2 2011/05/17 15:36:01 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:01 $
//
// Original author KLG based on Mu2eWorld constructPS
//
// Notes:
// Construct the PS. Parent volume is the air inside of the hall.

// C++ includes
#include <iostream>

// Mu2e includes.
#include "BeamlineGeom/inc/Beamline.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/constructPS.hh"
#include "Mu2eG4/inc/nestTubs.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"


using namespace std;

namespace mu2e {

  void constructPS( const VolumeInfo& parent, 
                    SimpleConfig const * const _config,
                    G4ThreeVector&  _primaryProtonGunOrigin,
                    G4RotationMatrix& _primaryProtonGunRotation
                    ){

    // Extract information from the config file.
    
    GeomHandle<Beamline> beamg;
    double solenoidOffset = beamg->solenoidOffset();
    double rTorus         = beamg->getTS().torusRadius();
    double ts1HalfLength  = beamg->getTS().getTS1().getHalfLength();

    double ps1HalfLength     = _config->getDouble("toyPS1.vacHalfLength");

    // Build the barrel of the cryostat.
    TubsParams psCryoParams( _config->getDouble("toyPS.rIn"),
                             _config->getDouble("toyPS.rOut"),
                             _config->getDouble("toyPS.CryoHalfLength"));

    MaterialFinder materialFinder(*_config);
    G4Material* psCryoMaterial = materialFinder.get("toyPS.materialName");

    // In the Mu2e coordinate system.
    double psCryoZ0 = -rTorus + -2.*ts1HalfLength - psCryoParams.zHalfLength;
    G4ThreeVector psCryoPosition( solenoidOffset, 0., psCryoZ0 );
    
    bool toyPSVisible        = _config->get<bool>("toyPS.visible",true);
    bool toyPSSolid          = _config->get<bool>("toyPS.solid",true);
    bool forceAuxEdgeVisible = _config->get<bool>("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config->get<bool>("g4.doSurfaceCheck",false);
    bool const placePV       = true;
    
    G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();

    // Toy model of the PS coils + cryostat. It needs real structure.
    VolumeInfo psCryoInfo = nestTubs( "PSCryo",
                                      psCryoParams,
                                      psCryoMaterial,
                                      0,
                                      psCryoPosition-_hallOriginInMu2e,
                                      parent,
                                      0,
                                      toyPSVisible,
                                      G4Color::Cyan(),
                                      toyPSSolid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

    // Build the main PS vacuum body.
    TubsParams ps1VacParams( 0.,
                             _config->getDouble("toyPS.rIn"),
                             ps1HalfLength);
    G4Material* vacuumMaterial  = materialFinder.get("toyDS.insideMaterialName");

    // Position in the Mu2e coordinate system.
    double ps1Z0     = -rTorus + -2.*ts1HalfLength - ps1HalfLength;
    G4ThreeVector ps1Position( solenoidOffset, 0., ps1Z0);    
    
    VolumeInfo ps1VacInfo   = nestTubs( "PS1Vacuum",
                                        ps1VacParams,
                                        vacuumMaterial,
                                        0,
                                        ps1Position-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        toyPSVisible,
                                        G4Colour::Green(),
                                        toyPSSolid,
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                        );
    
    // Build the production target.
    TubsParams prodTargetParams( 0.,
                                 _config->getDouble("targetPS_rOut"),
                                 _config->getDouble("targetPS_halfLength"));
    G4Material* prodTargetMaterial = materialFinder.get("targetPS_materialName");
    
    // Position in the Mu2e coordinate system.
    CLHEP::Hep3Vector prodTargetPosition = _config->getHep3Vector("productionTarget.position");

    // Rotation of production target.
    double targetPS_rotX = _config->getDouble("targetPS_rotX" )*CLHEP::degree;
    double targetPS_rotY = _config->getDouble("targetPS_rotY" )*CLHEP::degree;

    // G4 takes ownership of this G4RotationMatrix object.
    // Passive rotation. See Mu2e-doc-938.
    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
    G4RotationMatrix* prodTargetRotation = reg.add(G4RotationMatrix());
    prodTargetRotation->rotateY( -targetPS_rotY);
    prodTargetRotation->rotateX( -targetPS_rotX);

    bool prodTargetVisible = _config->get<bool>("targetPS.visible",true);
    bool prodTargetSolid   = _config->get<bool>("targetPS.solid",true);

    VolumeInfo prodTargetInfo   = nestTubs( "ProductionTarget",
                                            prodTargetParams,
                                            prodTargetMaterial,
                                            prodTargetRotation,
                                            prodTargetPosition-ps1Position,
                                            ps1VacInfo,
                                            0,
                                            prodTargetVisible,
                                            G4Colour::Magenta(),
                                            prodTargetSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );

    
    // Set the parameters of the transformation from the PrimaryProtonGun
    // coordinates to G4 coordinates.  This needs an active sense rotation,
    // the opposite of what G4 needed.
    _primaryProtonGunRotation = prodTargetRotation->inverse();

    G4ThreeVector prodTargetFaceLocal(0.,0.,prodTargetParams.zHalfLength);
    _primaryProtonGunOrigin = prodTargetPosition + VolumeInfo::getMu2eOriginInWorld() + 
      _primaryProtonGunRotation*prodTargetFaceLocal;
    
  } // end Mu2eWorld::constructPS
  }
