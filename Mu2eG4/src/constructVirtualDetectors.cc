//
// Free function to create the virtual detectors
//
// $Id: constructVirtualDetectors.cc,v 1.1 2011/01/05 21:04:47 genser Exp $
// $Author: genser $
// $Date: 2011/01/05 21:04:47 $
//
// Original author KLG based on Mu2eWorld constructVirtualDetectors
//
// C++ includes

#include <iostream>

// Mu2e includes.

#include "Mu2eG4/inc/constructVirtualDetectors.hh"

#include "BeamlineGeom/inc/Beamline.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/VirtualDetectorSD.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "TargetGeom/inc/Target.hh"
#include "VirtualDetectorGeom/inc/VirtualDetector.hh"

// G4 includes

#include "G4Material.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4Color.hh"


using namespace std;

namespace mu2e {

  // Construct the virtual detectors

  void constructVirtualDetectors( SimpleConfig const * const _config ){

    // Place virtual detectors

    bool vdVisible           = _config->getBool("vd.visible",true);
    bool vdSolid             = _config->getBool("vd.solid",true);
    bool forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
    bool const placePV       = true;
    
    GeomHandle<VirtualDetector> vdg;
    if( vdg->nDet()<=0 ) return;

    GeomHandle<Beamline> beamg;
    double rVac         = CLHEP::mm * beamg->getTS().innerRadius();

    GeomHandle<Target> target;
    double rTarget      = CLHEP::mm * target->cylinderRadius();

    double vdHalfLength = CLHEP::mm * vdg->getHalfLength();

    MaterialFinder materialFinder(*_config);
    G4Material* vacuumMaterial     = materialFinder.get("toyDS.insideMaterialName");

    TubsParams vdParams(0,rVac,vdHalfLength);
    TubsParams vdParamsTarget(0,rTarget,vdHalfLength);

    // Virtual Detectors 1 and 2 are placed inside TS1

    G4VSensitiveDetector* vdSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::VirtualDetector());

    G4Helper* _helper = &(*(edm::Service<G4Helper>()));

    for( int id=1; id<=2; ++id) if( vdg->exist(id) ) {
      VolumeInfo const & parent = _helper->locateVolInfo("ToyTS1Vacuum");
      ostringstream name;
      name << "VirtualDetector" << id;
      VolumeInfo vd = nestTubs( name.str(), vdParams, vacuumMaterial, 0,
                                vdg->getLocal(id),
                                parent,
                                id, vdVisible, G4Color::Red(), vdSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck );
      vd.logical->SetSensitiveDetector(vdSD);
    }

    // Virtual Detectors 3-6 are placed inside TS3

    for( int id=3; id<=6; ++id) if( vdg->exist(id) ) {
      VolumeInfo const & parent = _helper->locateVolInfo("ToyTS3Vacuum");
      ostringstream name;
      name << "VirtualDetector" << id;
      VolumeInfo vd = nestTubs( name.str(), vdParams, vacuumMaterial, 0,
                                vdg->getLocal(id),
                                parent,
                                id, vdVisible, G4Color::Red(), vdSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck );
      vd.logical->SetSensitiveDetector(vdSD);
    }

    // Virtual Detectors 7-8 are placed inside TS5

    for( int id=7; id<=8; ++id) if( vdg->exist(id) ) {
      VolumeInfo const & parent = _helper->locateVolInfo("ToyTS5Vacuum");
      ostringstream name;
      name << "VirtualDetector" << id;
      VolumeInfo vd = nestTubs( name.str(), vdParams, vacuumMaterial, 0,
                                vdg->getLocal(id),
                                parent,
                                id, vdVisible, G4Color::Red(), vdSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck );
      vd.logical->SetSensitiveDetector(vdSD);
    }

  }

}
