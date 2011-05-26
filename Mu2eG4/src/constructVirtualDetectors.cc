//
// Free function to create the virtual detectors
//
// $Id: constructVirtualDetectors.cc,v 1.6 2011/05/26 22:09:53 genser Exp $
// $Author: genser $
// $Date: 2011/05/26 22:09:53 $
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
#include "VirtualDetectorGeom/inc/VirtualDetector.hh"
#include "TTrackerGeom/inc/TTracker.hh"

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

    int static const verbosityLevel = _config->getInt("vd.verbosityLevel",0);

    bool vdVisible           = _config->getBool("vd.visible",true);
    bool vdSolid             = _config->getBool("vd.solid",true);
    bool forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
    bool const placePV       = true;

    GeomHandle<VirtualDetector> vdg;
    if( vdg->nDet()<=0 ) return;

    GeomHandle<Beamline> beamg;
    double rVac         = CLHEP::mm * beamg->getTS().innerRadius();

    double vdHalfLength = CLHEP::mm * vdg->getHalfLength();

    MaterialFinder materialFinder(*_config);
    G4Material* vacuumMaterial     = materialFinder.get("toyDS.insideMaterialName");

    TubsParams vdParams(0,rVac,vdHalfLength);

    // Virtual Detectors 1 and 2 are placed inside TS1

    G4VSensitiveDetector* vdSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::VirtualDetector());

    G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));

    for( int vdid=1; vdid<=2; ++vdid) if( vdg->exist(vdid) ) {
      VolumeInfo const & parent = _helper->locateVolInfo("ToyTS1Vacuum");
      ostringstream name;
      name << "VirtualDetector" << vdid;
      VolumeInfo vd = nestTubs( name.str(), vdParams, vacuumMaterial, 0,
                                vdg->getLocal(vdid),
                                parent,
                                vdid, vdVisible, G4Color::Red(), vdSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck );
      vd.logical->SetSensitiveDetector(vdSD);
    }

    // Virtual Detectors 3-6 are placed inside TS3

    for( int vdid=3; vdid<=6; ++vdid) if( vdg->exist(vdid) ) {
      VolumeInfo const & parent = _helper->locateVolInfo("ToyTS3Vacuum");
      ostringstream name;
      name << "VirtualDetector" << vdid;
      VolumeInfo vd = nestTubs( name.str(), vdParams, vacuumMaterial, 0,
                                vdg->getLocal(vdid),
                                parent,
                                vdid, vdVisible, G4Color::Red(), vdSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck );
      vd.logical->SetSensitiveDetector(vdSD);
    }

    // Virtual Detectors 7-8 are placed inside TS5

    for( int vdid=7; vdid<=8; ++vdid) if( vdg->exist(vdid) ) {
      VolumeInfo const & parent = _helper->locateVolInfo("ToyTS5Vacuum");
      ostringstream name;
      name << "VirtualDetector" << vdid;
      VolumeInfo vd = nestTubs( name.str(), vdParams, vacuumMaterial, 0,
                                vdg->getLocal(vdid),
                                parent,
                                vdid, vdVisible, G4Color::Red(), vdSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck );
      vd.logical->SetSensitiveDetector(vdSD);
    }

    // Virtual Detectors 9-10 are placed inside DS2, just before and after stopping target

    // If there is no neutron absorber, virtual detectors 9 and 10 extend to
    // inner wall of DS2 minus 5 mm. If neutron absorber is define, these
    // detectors extend to neutron absorber minus 5 mm.

    double Ravr = _config->getDouble("toyDS.rIn");
    double deltaR = 0;
    double Z0 = 0;
    double deltaZ = 1.0;

    if ( _config->getBool("hasNeutronAbsorber",false) ) {
      double NAIInnerRadius0     = _config->getDouble("neutronabsorber.internalInnerRadius0");
      double NAIInnerRadius1     = _config->getDouble("neutronabsorber.internalInnerRadius1");
      Ravr   = (NAIInnerRadius0+NAIInnerRadius1)/2;
      deltaR = (NAIInnerRadius1-NAIInnerRadius0);
      Z0     = _config->getDouble("neutronabsorber.internalZ01");
      deltaZ = 2.0 *_config->getDouble("neutronabsorber.internalHalfLengthZ01");
    }

    for( int vdid=9; vdid<=10; ++vdid) if( vdg->exist(vdid) ) {
      
      double zvd = vdg->getGlobal(vdid).z();
      double rvd = Ravr + deltaR/deltaZ*(zvd-Z0) - 5.0;

      if ( verbosityLevel > 0) {
        cout << "VD " << vdid << " : " << zvd << " " << rvd << endl;
      }

      TubsParams vdParamsTarget(0.,rvd,vdHalfLength);

      VolumeInfo const & parent = _helper->locateVolInfo("ToyDS2Vacuum");
      ostringstream name;
      name << "VirtualDetector" << vdid;
      VolumeInfo vd = nestTubs( name.str(), vdParamsTarget, vacuumMaterial, 0,
                                vdg->getLocal(vdid),
                                parent,
                                vdid, vdVisible, G4Color::Red(), vdSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck );
      vd.logical->SetSensitiveDetector(vdSD);
    }

    // placing virtual detectors in the middle of the ttracker

    // check if ttracker exists and if the number of devices
    // ttracker.numDevices is even is done in VirtualDetectorMaker

    int vdid = 11;
    if( vdg->exist(vdid) ) {
      
      // the radius of tracker mother
      TTracker const & ttracker = *(GeomHandle<TTracker>());
      double orvd = ttracker.getTrackerEnvelopeParams().outerRadius();
      double irvd = ttracker.getTrackerEnvelopeParams().innerRadius();

      if ( verbosityLevel > 0) {
        double zvd = vdg->getGlobal(vdid).z();
        cout << "VD " << vdid << " z, r : " << zvd << " " << irvd << " " << orvd << endl;
      }

      TubsParams vdParamsTTracker(irvd,orvd,vdHalfLength);

      VolumeInfo const & parent = _helper->locateVolInfo("TrackerMother");
      ostringstream name;
      name << "VirtualDetector" << vdid;
      VolumeInfo vd = nestTubs( name.str(), vdParamsTTracker, vacuumMaterial, 0,
                                vdg->getLocal(vdid),
                                parent,
                                vdid, vdVisible, G4Color::Red(), vdSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck );
      vd.logical->SetSensitiveDetector(vdSD);

      vdid = 12;
      if( vdg->exist(vdid) ) {

        TubsParams vdParamsTTrackerInner(0.,irvd,vdHalfLength);

      // VD 12 is placed inside the ttracker at the same z position as
      // VD 11 but from radius 0 to the inner radius of the ttracker
      // mother volume. However, its mother volume is ToyDS3Vacuum
      // which has a different offset. We will use the global offset
      // here (!) as DS is not in the geometry service yet

        VolumeInfo const & parent = _helper->locateVolInfo("ToyDS3Vacuum");

        G4ThreeVector vdLocalOffset = vdg->getGlobal(vdid) - parent.centerInMu2e();

        if ( verbosityLevel > 0) {
          double zvd = vdg->getGlobal(vdid).z();
          cout << "VD " << vdid << " z, r : " << zvd  << " " << irvd << endl;
        }

        ostringstream name;
        name << "VirtualDetector" << vdid;
        VolumeInfo vd = nestTubs( name.str(), vdParamsTTrackerInner, vacuumMaterial, 0,
                                  vdLocalOffset,
                                  parent,
                                  vdid, vdVisible, G4Color::Red(), vdSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck );
        vd.logical->SetSensitiveDetector(vdSD);
      }

    }

  }

}
