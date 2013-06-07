//
// Free function to create Neutron Absorbers in G4
//
// $Id: constructInternalNeutronAbsorber.cc,v 1.4 2013/06/07 17:43:30 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/06/07 17:43:30 $
//
// Original author KLG
//
// Notes:
// Construct the InternalNeutron Absorbers in G4


// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Mu2e includes.

#include "Mu2eG4/inc/constructInternalNeutronAbsorber.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "InternalNeutronAbsorberGeom/inc/InternalNeutronAbsorber.hh"
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

  void constructInternalNeutronAbsorber(SimpleConfig const & _config){

    int const verbosityLevel = _config.getInt("neutronabsorber.verbosityLevel",0);

    // Extract information from the config file.
    bool NAVisible             = _config.getBool("neutronabsorber.visible");
    bool NASolid               = _config.getBool("neutronabsorber.solid");

    bool const forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    // Access to the G4HelperService.
    G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));

    // now constructing the internal neutron absorber
    // it is placed inside DS2Vacuum & DS3Vacuum like the protonabs1 & 2

    // Fetch InternalNeutronAborber geometry
    GeomHandle<InternalNeutronAbsorber> nai;

    double const NAI2FrontZ = nai->position().z() - nai->halfLength2();
    double const NAI1FrontZ = nai->position().z() - nai->halfLength1();

    if ( verbosityLevel > 0) {
      cout << __func__ << " NAI1 FrontZ                       : " << NAI1FrontZ    << endl;
      cout << __func__ << " NAI2 FrontZ                       : " << NAI2FrontZ    << endl;
      cout << __func__ << " NAI1 ExtentZ                      : " << NAI2FrontZ + 2*nai->halfLength2() << endl;
      cout << __func__ << " NAI1 ExtentZ                      : " << NAI1FrontZ + 2*nai->halfLength1() << endl;
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

    if ( ds2FrontZ>NAI1FrontZ ){
      throw cet::exception("GEOM")
        << "Internal Neutron Absorber 1 front is outside of DS2Vaccum \n";
    }

    if ( ds->vac_zLocDs23Split()<NAI2FrontZ+2*nai->halfLength2() ||
         ds->vac_zLocDs23Split()<NAI1FrontZ+2*nai->halfLength1() ){
      throw cet::exception("GEOM")
        << "Internal Neutron Absorber is outside of DS2Vaccum \n";
    }

    // Get DS2Vacuum & DS3Vacuum info
    VolumeInfo const & ds2VacuumInfo = _helper->locateVolInfo("DS2Vacuum");

    CLHEP::Hep3Vector NAIOffset  =  nai->position() - ds2VacuumInfo.centerInMu2e();

    if ( verbosityLevel > 0) {
      cout << __func__ << " DS2VacuumInfo.centerInMu2e()  : " << ds2VacuumInfo.centerInMu2e() << endl;
      cout << __func__ << " NAIOffset                       : "  << NAIOffset  << endl;
    }

    TubsParams NAI2Params( nai->rIn2(), nai->rOut(), nai->halfLength2() );
    TubsParams NAI1Params( nai->rIn1(), nai->rIn2(), nai->halfLength1() );

    G4Material* NAI2Material = findMaterialOrThrow( nai->material2() );
    G4Material* NAI1Material = findMaterialOrThrow( nai->material1() );

    nestTubs("InternalNeutronAbsorber2",
             NAI2Params,
             NAI2Material,
             0,
             NAIOffset,
             ds2VacuumInfo,
             0,
             NAVisible,
             G4Colour::Cyan(),
             NASolid,
             forceAuxEdgeVisible,
             placePV,
             doSurfaceCheck
             );

    nestTubs("InternalNeutronAbsorber1",
             NAI1Params,
             NAI1Material,
             0,
             NAIOffset,
             ds2VacuumInfo,
             0,
             NAVisible,
             G4Colour::Cyan(),
             NASolid,
             forceAuxEdgeVisible,
             placePV,
             doSurfaceCheck
             );

  } // end of constructInternalNeutronAbsorber;

}
