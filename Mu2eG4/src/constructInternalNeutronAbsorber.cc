//
// Free function to create Neutron Absorbers in G4
//
// $Id: constructInternalNeutronAbsorber.cc,v 1.1 2013/03/28 13:02:38 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/03/28 13:02:38 $
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

    // the internal neutron absorber is placed in the DS2/3 vacuum
    // split accordingly to deal with the ds2/ds3 split

    // Extract information from the config file.
    bool NAVisible             = _config.getBool("neutronabsorber.visible");
    bool NASolid               = _config.getBool("neutronabsorber.solid");

    bool const forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    // Access to the G4HelperService.
    G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));

    // now constructing the internal neutron absorber
    // it is placed inside ToyDS2Vacuum & ToyDS3Vacuum like the protonabs1 & 2

    // Fetch InternalNeutronAborber geometry
    GeomHandle<InternalNeutronAbsorber> nai;

    double const NAI1FrontZ     = nai->positionAbs1().z() - nai->halfLengthAbs1();
    double const NAI12BoundaryZ = nai->positionAbs1().z() + nai->halfLengthAbs1();

    if ( verbosityLevel > 0) {
      cout << __func__ << " NAI1OffsetInMu2e                 : " << nai->positionAbs1() << endl;
    }

    // halflength of second neutron absorber constrained by z-extent
    // of ToyDS3Vacuum

    // Fetch DS geom object
    GeomHandle<DetectorSolenoid> ds;

    double const ds2HalfLength = ds->halfLengthDs2();
    double const ds3FullLength = ds->coilZMax() - ds->zLocDs23Split();

    double const NAIHalfLengthZ02 = (ds->coilZMax() - NAI12BoundaryZ)*0.5;

    double const NAIZ02    = nai->positionAbs1().z() + nai->halfLengthAbs1() + NAIHalfLengthZ02;
    double const NAI2EndZ  = NAIZ02 + NAIHalfLengthZ02;

    if ( verbosityLevel > 0) {
      cout << __func__ << " NAIZ01                           : " << nai->positionAbs1().z()  << endl;
      cout << __func__ << " NAIZ02                           : " << NAIZ02  << endl;
      cout << __func__ << " NAI1FrontZ                       : " << NAI1FrontZ    << endl;
      cout << __func__ << " NAI12BoundaryZ                   : " << NAI12BoundaryZ << endl;
      cout << __func__ << " NAI2EndZ                         : " << NAI2EndZ      << endl;
    }

    // the Z offset of the InternalNeutronAbsorber2

    // we need to calculate where is the boundary of
    // InternalNeutronAbsorber1,2 wrt boudary between ToyDS2Vacuum &
    // ToyDS3Vacuum

    if ( verbosityLevel > 0) {
      cout << __func__ << " toyDS2Vacuum half-length                 : " <<  ds2HalfLength  << endl;
      cout << __func__ << " toyDS3Vacuum full length in DS (polycone): " <<  ds3FullLength  << endl;
    }

    // what are the offsets of the DS2,3Vacuum in the Mu2e system?
    // the ds2/ds3 boundary in mu2e coordinates

    double const ds2FrontZG = ds->zLocDs23Split() - 2.*ds->halfLengthDs2();

    if ( verbosityLevel > 0) {
      cout << __func__ << " ds2FrontZG                       : " << ds2FrontZG    << endl;
      cout << __func__ << " ds23BoundaryG                    : " << ds->zLocDs23Split() << endl;
      cout << __func__ << " ds3EndZG                         : " << ds->coilZMax()      << endl;
    }

    // certain combinations are illegal; we will assume that the
    // conical absorber is fully contained in ds2

    if ( ds2FrontZG>NAI1FrontZ ){
      throw cet::exception("GEOM")
        << "Internal Neutron Absorber 1 front is outside of DS2Vaccum \n";
    }

    if ( ds->zLocDs23Split()<NAI12BoundaryZ ){
      throw cet::exception("GEOM")
        << "Internal Neutron Absorber 1 end is outside of DS2Vaccum \n";
    }

    if ( ds->coilZMax()<NAI2EndZ ){
      throw cet::exception("GEOM")
        << "Internal Neutron Absorber 2 end is outside of DS3Vaccum \n";
    }

    // tube parameters for InternalNeutronAbsorber2

    // we need to construct two InternalNeutronAbsorber2 one in ds2 and another one in ds3
    // we need to calculate their offsets
    // x-offsets are inline with those of InternalNeutronAbsorber1

    double const halfDeltaBoundary = (ds->zLocDs23Split() - NAI12BoundaryZ)*.5;
    double const averageBoundary   = (ds->zLocDs23Split() + NAI12BoundaryZ)*.5;

    CLHEP::Hep3Vector NAI23OffsetInMu2e =
      CLHEP::Hep3Vector(nai->positionAbs1().x(),0.,NAIZ02+halfDeltaBoundary);

    if ( verbosityLevel > 0) {
      cout << __func__ << " NAI23OffsetInMu2e                : " << NAI23OffsetInMu2e << endl;
    }

    CLHEP::Hep3Vector NAI22OffsetInMu2e =
      CLHEP::Hep3Vector(nai->positionAbs1().x(),0.,averageBoundary);

    if ( verbosityLevel > 0) {
      cout << __func__ << " NAI22OffsetInMu2e                : " << NAI22OffsetInMu2e << endl;
    }

    // Get ToyDS2Vacuum & ToyDS3Vacuum info
    VolumeInfo const & toyDS2VacuumInfo = _helper->locateVolInfo("ToyDS2Vacuum");
    VolumeInfo const & toyDS3VacuumInfo = _helper->locateVolInfo("ToyDS3Vacuum");

    CLHEP::Hep3Vector NAI1Offset  =  nai->positionAbs1() - toyDS2VacuumInfo.centerInMu2e();

    // constructing InternalNeutronAbsorber1 as concentric annular cylinders
    vector<TubsParams> NAI1Params;
    vector<VolumeInfo> NAI1VolumeInfo;
    for ( size_t i = 0 ; i < nai->rInAbs1Vec().size() ; i++ ) {
      if ( i < nai->rInAbs1Vec().size()-1 ) // do inner cylinders
        NAI1Params.push_back( TubsParams( nai->rInAbs1Vec().at(i), nai->rInAbs1Vec().at(i+1), nai->halfLengthAbs1()) );
      else                           // do outer cylinder
        NAI1Params.push_back( TubsParams( nai->rInAbs1Vec().back(), nai->rOut(), nai->halfLengthAbs1() ) );
        
      ostringstream detname( "InternalNeutronAbsorber1_");
      detname << i;
	  
      G4Material* NAI1Material = findMaterialOrThrow( nai->materialAbs1Vec().at(i) );
      NAI1VolumeInfo.push_back( nestTubs(detname.str(),
                                         NAI1Params.back(),
                                         NAI1Material,
                                         0,
                                         NAI1Offset,
                                         toyDS2VacuumInfo,
                                         0,
                                         NAVisible,
                                         G4Colour::Cyan(),
                                         NASolid,
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         ) );
	  
    }

    //constructing split internal neutron absorber 2
    CLHEP::Hep3Vector NAI22Offset =  NAI22OffsetInMu2e - toyDS2VacuumInfo.centerInMu2e();
    CLHEP::Hep3Vector NAI23Offset =  NAI23OffsetInMu2e - toyDS3VacuumInfo.centerInMu2e();

    if ( verbosityLevel > 0) {
      cout << __func__ << " toyDS2VacuumInfo.centerInMu2e()  : " << toyDS2VacuumInfo.centerInMu2e() << endl;
      cout << __func__ << " toyDS3VacuumInfo.centerInMu2e()  : " << toyDS3VacuumInfo.centerInMu2e() << endl;
      cout << __func__ << " NAI1Offset                       : " << NAI1Offset  << endl;
      cout << __func__ << " NAI22Offset                      : " << NAI22Offset << endl;
      cout << __func__ << " NAI23Offset                      : " << NAI23Offset << endl;
    }

    TubsParams NAI22Params( nai->rInAbs2(),
                            nai->rOut(),
                            halfDeltaBoundary);

    TubsParams NAI23Params( nai->rInAbs2(),
                            nai->rOut(),
                            NAIHalfLengthZ02-halfDeltaBoundary);


    G4Material* NAMaterial = findMaterialOrThrow( nai->materialAbs2() );

    VolumeInfo INA22Info  = nestTubs("InternalNeutronAbsorber22",
                                     NAI22Params,
                                     NAMaterial,
                                     0,
                                     NAI22Offset,
                                     toyDS2VacuumInfo,
                                     0,
                                     NAVisible,
                                     G4Colour::Cyan(),
                                     NASolid,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     doSurfaceCheck
                                     );

    VolumeInfo INA23Info  = nestTubs("InternalNeutronAbsorber23",
                                     NAI23Params,
                                     NAMaterial,
                                     0,
                                     NAI23Offset,
                                     toyDS3VacuumInfo,
                                     0,
                                     NAVisible,
                                     G4Colour::Cyan(),
                                     NASolid,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     doSurfaceCheck
                                     );

    if ( verbosityLevel > 0) {
      double zhl         = static_cast<G4Tubs*>(INA23Info.solid)->GetZHalfLength();
      CLHEP::Hep3Vector const & IntNeutronAbs23OffsetInMu2e = INA23Info.centerInMu2e();
      double IntNeutronAbs23OffsetInMu2eZ = IntNeutronAbs23OffsetInMu2e[CLHEP::Hep3Vector::Z];
      cout << __func__ << " INA23Info Z center in Mu2e    : " << IntNeutronAbs23OffsetInMu2eZ << endl;
      cout << __func__ << " INA23Info Z extent in Mu2e    : " <<
        IntNeutronAbs23OffsetInMu2eZ - zhl << ", " << IntNeutronAbs23OffsetInMu2eZ + zhl << endl;
    } 

  } // end of constructNeutronAbsorber;

}
