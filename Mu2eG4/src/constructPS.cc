//
// Free function to create  Production Solenoid and Production Target.
//
// $Id: constructPS.cc,v 1.8 2012/02/27 06:05:35 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/27 06:05:35 $
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
#include "ProductionTargetGeom/inc/ProductionTarget.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/constructPS.hh"
#include "Mu2eG4/inc/nestTubs.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"

#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

namespace mu2e {

  void constructPS(const VolumeInfo& parent, const SimpleConfig * const _config) {
    
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
    double psCryoZ0 = -rTorus + -2.*ts1HalfLength - psCryoParams.zHalfLength();
    G4ThreeVector psCryoPosition( solenoidOffset, 0., psCryoZ0 );

    bool toyPSVisible        = _config->getBool("toyPS.visible",true);
    bool toyPSSolid          = _config->getBool("toyPS.solid",true);
    bool forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
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
    GeomHandle<ProductionTarget> tgt;
    TubsParams prodTargetParams( 0., tgt->rOut(), tgt->halfLength());

    G4Material* prodTargetMaterial = materialFinder.get("targetPS_materialName");
    bool prodTargetVisible = _config->getBool("targetPS.visible",true);
    bool prodTargetSolid   = _config->getBool("targetPS.solid",true);

    VolumeInfo prodTargetInfo   = nestTubs( "ProductionTarget",
                                            prodTargetParams,
                                            prodTargetMaterial,
                                            &tgt->productionTargetRotation(),
                                            tgt->position() - ps1Position,
                                            ps1VacInfo,
                                            0,
                                            prodTargetVisible,
                                            G4Colour::Magenta(),
                                            prodTargetSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );


    // To compare with g4beamline studies: close the vacuum with a solid disk
    const double toyEnclosureThickness(_config->getDouble("toyPS.toyEnclosure.Thickness", 0)*CLHEP::mm);
    if(toyEnclosureThickness > 0) {
      
      TubsParams diskParam(0, _config->getDouble("toyPS.rOut"), 0.5*toyEnclosureThickness);
      CLHEP::Hep3Vector toyEnclosurePosition = psCryoPosition 
	- CLHEP::Hep3Vector(0, 0, psCryoParams.zHalfLength())
	- CLHEP::Hep3Vector(0, 0, 0.5*toyEnclosureThickness)
	;

      nestTubs( "PS1ToyEnclosure",
		diskParam,
		materialFinder.get("toyPS.toyEnclosure.materialName"),
		0,
		toyEnclosurePosition - _hallOriginInMu2e,
		parent,
		0,
		_config->getBool("toyPS.toyEnclosure.visible", true),
		G4Colour::Magenta(),
		_config->getBool("toyPS.toyEnclosure.solid", true),
		forceAuxEdgeVisible,
		placePV,
		doSurfaceCheck
		);
    }

  } // end Mu2eWorld::constructPS
}
