//
// Free function to create TSdA Neutron Absorber in G4
//
// $Id: constructTSdA.cc,v 1.6 2013/10/12 00:19:47 brownd Exp $
// $Author: brownd $
// $Date: 2013/10/12 00:19:47 $
//
// Original author KLG
//
// Notes:
// Construct the Internal Neutron Absorbers in G4
// David Norvil Brown (the other one):  rename to TSdA for consistency
// with TDR, and update - May 2015.

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Mu2e includes.

#include "Mu2eG4/inc/constructTSdA.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/StraightSection.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "BeamlineGeom/inc/TSdA.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4VSolid.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4VPhysicalVolume.hh"

using namespace std;

namespace mu2e {

  void constructTSdA(SimpleConfig const & _config){

    int const verbosityLevel = _config.getInt("tsda.verbosityLevel",0);
    double tmpRin = _config.getDouble("tsda.rin",0.0)*CLHEP::mm;


    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "tsda", "tsda");

    bool const NAVisible           = geomOptions->isVisible("tsda"); 
    bool const NASolid             = geomOptions->isSolid("tsda"); 
    const bool doSurfaceCheck      = geomOptions->doSurfaceCheck("tsda"); 
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible("tsda"); 
    const bool placePV             = geomOptions->placePV("tsda");  
 
    // Access to the G4HelperService.
    G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));

    // now constructing the internal neutron absorber
    // it is placed inside DS2Vacuum & DS3Vacuum like the protonabs1 & 2

    // Fetch TSdA geometry
    GeomHandle<TSdA> atsd;

    if ( tmpRin < 1.0e-06 ) {  // just a check for zero
      // TS5 outer radius 
      GeomHandle<Beamline> beamg;
      const StraightSection * ts5out = beamg->getTS().getTSCryo<StraightSection>( TransportSolenoid::TSRegion::TS5,
										  TransportSolenoid::TSRadialPart::IN );
      tmpRin = ts5out->rIn();
    }
    if ( verbosityLevel > 0) {
      cout << __func__ << " TSdA rin                      : "<< tmpRin << endl;
      cout << __func__ << " TSdA rout			  : "<< atsd->r4() << endl;
    }

    // we need to calculate where the DS2Vacuum volume is;
    // the int. neutron absorber must fit inside it

    // Fetch DS geom object
    GeomHandle<DetectorSolenoid> ds;

    double const ds2FrontZ = ds->vac_zLocDs23Split() - 2.*ds->vac_halfLengthDs2();
    if ( verbosityLevel > 0) {
      cout << __func__ << " DS2Vacuum extent               : [" 
           << ds2FrontZ << " , " 
           << ds->vac_zLocDs23Split() << endl;
    }

    // certain combinations are illegal; we will assume that the
    // conical absorber is fully contained in ds2


    // Get DS2Vacuum & DS3Vacuum info
    VolumeInfo const & ds2VacuumInfo = _helper->locateVolInfo("DS2Vacuum");
    VolumeInfo const & ds3VacuumInfo = _helper->locateVolInfo("DS3Vacuum");

    CLHEP::Hep3Vector ATSDOffset = atsd->position() - ds2VacuumInfo.centerInMu2e();

    CLHEP::Hep3Vector ATSD4Offset(0.0,0.0,-ds->vac_halfLengthDs2() + atsd->halfLength4());

    if ( verbosityLevel > 0) {
      cout << __func__ << " DS2VacuumInfo.centerInMu2e()  : " << ds2VacuumInfo.centerInMu2e() << endl;
      cout << __func__ << " DS3VacuumInfo.centerInMu2e()  : " << ds3VacuumInfo.centerInMu2e() << endl;
      cout << __func__ << " TSdA Offset                     : "  << ATSDOffset  << endl;
      cout << __func__ << " TSdA Offset2			  : "  << ATSD4Offset  << endl;
      cout << __func__ << " TSdA center in Mu2e		  : "  << (ds2VacuumInfo.centerInMu2e() + ATSD4Offset).z()  << endl;
    }

    TubsParams ATSD4Params( tmpRin, atsd->r4(), atsd->halfLength4() );

    G4Material* ATSD4Material = findMaterialOrThrow( atsd->material4() );

    nestTubs("TSdA4",
             ATSD4Params,
             ATSD4Material,
             0,
             ATSD4Offset,
             ds2VacuumInfo,
             0,
             NAVisible,
             G4Colour::Cyan(),
             NASolid,
             forceAuxEdgeVisible,
             placePV,
             doSurfaceCheck
             );


  } // end of constructTSdA;

}
