//
// Free function to create Neutron Absorbers in G4
//
// $Id: constructInternalNeutronAbsorber.cc,v 1.6 2013/10/12 00:19:47 brownd Exp $
// $Author: brownd $
// $Date: 2013/10/12 00:19:47 $
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
#include "BeamlineGeom/inc/StraightSection.hh"
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

    int const verbosityLevel = _config.getInt("intneutronabs.verbosityLevel",0);

    // Extract information from the config file.
    bool NAVisible             = _config.getBool("intneutronabs.visible");
    bool NASolid               = _config.getBool("intneutronabs.solid");

    bool const forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    // Access to the G4HelperService.
    G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));

    // now constructing the internal neutron absorber
    // it is placed inside DS2Vacuum & DS3Vacuum like the protonabs1 & 2

    // Fetch InternalNeutronAborber geometry
    GeomHandle<InternalNeutronAbsorber> nai;

    double const NAI3FrontZ = nai->position().z() + nai->halfLength2();
    double const NAI2FrontZ = nai->position().z() - nai->halfLength2();
    double const NAI1FrontZ = nai->position().z() - nai->halfLength1();

  // TS5 outer radius 
    GeomHandle<Beamline> beamg;
    const StraightSection * ts5out = beamg->getTS().getTSCryo<StraightSection>( TransportSolenoid::TSRegion::TS5,
                                                                                TransportSolenoid::TSRadialPart::IN );

    if ( verbosityLevel > 0) {
      cout << __func__ << " NAI1 FrontZ                       : " << NAI1FrontZ    << endl;
      cout << __func__ << " NAI2 FrontZ                       : " << NAI2FrontZ    << endl;
      cout << __func__ << " NAI3 FrontZ                       : " << NAI3FrontZ    << endl;
      cout << __func__ << " NAI1 ExtentZ                      : " << NAI1FrontZ + 2*nai->halfLength1() << endl;
      cout << __func__ << " NAI2 ExtentZ                      : " << NAI2FrontZ + 2*nai->halfLength2() << endl;
      cout << __func__ << " NAI3 ExtentZ                      : " << NAI3FrontZ + 2*nai->halfLength3() << endl;
      cout << __func__ << " NAI3 rin1                      : " << nai->rIn1() << endl;
      cout << __func__ << " NAI3 rin2                      : " << nai->rIn2() << endl;
      cout << __func__ << " NAI3 rin3                      : " << nai->rIn3() << endl;
      cout << __func__ << " NAI4 rin                      : " << ts5out->rIn() << endl;
      cout << __func__ << " NAI4 rout			  : " << nai->r4() << endl;
      cout << __func__ << " NAI3 rout                      : " << nai->rOut() << endl;
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
    VolumeInfo const & ds3VacuumInfo = _helper->locateVolInfo("DS3Vacuum");

    // downstream section is displaced, and split across DS2 and DS3
    double NAI3zedges[2] = {nai->position().z() + nai->halfLength2(),
			nai->position().z() + nai->halfLength2()+ 2*nai->halfLength3()};
    double NAI3azcenter = 0.5*(NAI3zedges[0] + ds->vac_zLocDs23Split());
    double NAI3ahalflength = 0.5*(ds->vac_zLocDs23Split()-NAI3zedges[0]);

    double NAI3bzcenter = 0.5*(ds->vac_zLocDs23Split()+NAI3zedges[1]);
    double NAI3bhalflength = 0.5*(NAI3zedges[1] - ds->vac_zLocDs23Split());

    CLHEP::Hep3Vector NAIOffset = nai->position() - ds2VacuumInfo.centerInMu2e();

    CLHEP::Hep3Vector NAI3acenter = nai->position();
    NAI3acenter.setZ(NAI3azcenter);
    CLHEP::Hep3Vector NAI3aOffset  = NAI3acenter - ds2VacuumInfo.centerInMu2e();

    CLHEP::Hep3Vector NAI3bcenter = nai->position();
    NAI3bcenter.setZ(NAI3bzcenter);
    CLHEP::Hep3Vector NAI3bOffset  = NAI3bcenter - ds3VacuumInfo.centerInMu2e();

    CLHEP::Hep3Vector NAI4Offset(0.0,0.0,-ds->vac_halfLengthDs2() + nai->halfLength4());

    if ( verbosityLevel > 0) {
      cout << __func__ << " DS2VacuumInfo.centerInMu2e()  : " << ds2VacuumInfo.centerInMu2e() << endl;
      cout << __func__ << " DS3VacuumInfo.centerInMu2e()  : " << ds3VacuumInfo.centerInMu2e() << endl;
      cout << __func__ << " NAIOffset                     : "  << NAIOffset  << endl;
      cout << __func__ << " NAI3aOffset                   : "  << NAI3aOffset  << endl;
      cout << __func__ << " NAI3bOffset                   : "  << NAI3bOffset  << endl;
      cout << __func__ << " NAI4Offset			  : "  << NAI4Offset  << endl;
      cout << __func__ << " NAI3a half-length             : "  << NAI3ahalflength  << endl;
      cout << __func__ << " NAI3b half-length             : "  << NAI3bhalflength  << endl;
    }

    TubsParams NAI4Params( ts5out->rIn(), nai->r4(), nai->halfLength4() );
    TubsParams NAI3aParams( nai->rIn3(), nai->rOut(), NAI3ahalflength );
    TubsParams NAI3bParams( nai->rIn3(), nai->rOut(), NAI3bhalflength );
    TubsParams NAI2Params( nai->rIn2(), nai->rOut(), nai->halfLength2() );
    TubsParams NAI1Params( nai->rIn1(), nai->rIn2(), nai->halfLength1() );

    G4Material* NAI4Material = findMaterialOrThrow( nai->material4() );
    G4Material* NAI3Material = findMaterialOrThrow( nai->material3() );
    G4Material* NAI2Material = findMaterialOrThrow( nai->material2() );
    G4Material* NAI1Material = findMaterialOrThrow( nai->material1() );

    nestTubs("InternalNeutronAbsorber4",
             NAI4Params,
             NAI4Material,
             0,
             NAI4Offset,
             ds2VacuumInfo,
             0,
             NAVisible,
             G4Colour::Cyan(),
             NASolid,
             forceAuxEdgeVisible,
             placePV,
             doSurfaceCheck
             );

    nestTubs("InternalNeutronAbsorber3a",
             NAI3aParams,
             NAI3Material,
             0,
             NAI3aOffset,
             ds2VacuumInfo,
             0,
             NAVisible,
             G4Colour::Cyan(),
             NASolid,
             forceAuxEdgeVisible,
             placePV,
             doSurfaceCheck
             );

    nestTubs("InternalNeutronAbsorber3b",
             NAI3bParams,
             NAI3Material,
             0,
             NAI3bOffset,
             ds3VacuumInfo,
             0,
             NAVisible,
             G4Colour::Cyan(),
             NASolid,
             forceAuxEdgeVisible,
             placePV,
             doSurfaceCheck
             );

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
