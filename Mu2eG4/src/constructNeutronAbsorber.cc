//
// Free function to create Neutron Absorbers in G4
//
// $Id: constructNeutronAbsorber.cc,v 1.2 2011/04/29 17:44:15 genser Exp $
// $Author: genser $
// $Date: 2011/04/29 17:44:15 $
//
// Original author KLG 
//
// Notes:
// Construct the Neutron Absorbers in G4


// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Mu2e includes.

#include "Mu2eG4/inc/constructNeutronAbsorber.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "G4Helper/inc/G4Helper.hh"
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
#include "G4SubtractionSolid.hh"
#include "G4VPhysicalVolume.hh"

using namespace std;

namespace mu2e {

  void constructNeutronAbsorber(SimpleConfig const * const _config){

    int const verbosityLevel = _config->getInt("neutronabsorber.verbosityLevel",0);

    // the absorber is split into two major pieces internal & external (wrt to the coil)

    // the internal part is placed in the DS2/3 vacuum
    // split accordingly to deal with the ds2/ds3 split

    // the external part is placed in the hall air (and does not need to be split)

    // Extract information from the config file.

    string NAMaterialName      = _config->getString("neutronabsorber.materialName");
    // string NADopantElementName = _config->getString("neutronabsorber.dopantElementName");
    // double NAMaterialDensity   = _config->getDouble("neutronabsorber.materialDensity");
    // double NADopantFraction    = _config->getDouble("neutronabsorber.dopantFraction");
    double NAIOuterRadius      = _config->getDouble("neutronabsorber.internalOuterRadius");
    // increasing inner radius
    double NAIInnerRadius0     = _config->getDouble("neutronabsorber.internalInnerRadius0");
    double NAIInnerRadius1     = _config->getDouble("neutronabsorber.internalInnerRadius1");
    double NAIInnerRadius2     = _config->getDouble("neutronabsorber.internalInnerRadius2");
    double NAIHalfLengthZ01    = _config->getDouble("neutronabsorber.internalHalfLengthZ01");
    double NAIHalfLengthZ02    = _config->getDouble("neutronabsorber.internalHalfLengthZ02");
    double NAIZ01              = _config->getDouble("neutronabsorber.internalZ01");
    double NAEHalfLengthZ      = _config->getDouble("neutronabsorber.externalHalfLengthZ");
    double NAEHalfLengthXY     = _config->getDouble("neutronabsorber.externalHalfLengthXY");
    double NAEHalfThickness    = _config->getDouble("neutronabsorber.externalHalfThickness");
    double NAEZ0               = _config->getDouble("neutronabsorber.externalZ0");

    bool NAVisible             = _config->getBool("neutronabsorber.visible");
    bool NASolid               = _config->getBool("neutronabsorber.solid");

    bool const forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    G4Material* NAMaterial = findMaterialOrThrow(NAMaterialName);

    // constructing External Neutron Absorber, placing it in HallAir (the name is hardcoded here...)

    // Access to the G4HelperService.
    G4Helper* _helper = &(*(edm::Service<G4Helper>()));
    
    VolumeInfo const & hallInfo = _helper->locateVolInfo("HallAir");

    // NAEZ0 is the offset in mu2e
    // the Z offset in World is NAEZ0+mu2eOrigin

    GeomHandle<Beamline> beamg;
    double solenoidOffset = -beamg->solenoidOffset(); 
    // this is an offset in X (and in what?) and should be negative?

    if ( verbosityLevel > 0) {
      cout << __func__ << " solenoidOffset                   : " << solenoidOffset  << endl;
    }

    CLHEP::Hep3Vector NAEOffsetInMu2e  = CLHEP::Hep3Vector(solenoidOffset,0.,NAEZ0);

    if ( verbosityLevel > 0) {
      cout << __func__ << " NAEOffsetInMu2e                  : " << NAEOffsetInMu2e << endl;
    }

    // now local offset in hallAir
    CLHEP::Hep3Vector NAEOffset =  NAEOffsetInMu2e - hallInfo.centerInMu2e();

    if ( verbosityLevel > 0) {
      cout << __func__ << " hallInfo.centerInMu2e()          : " << hallInfo.centerInMu2e() << endl;
      cout << __func__ << " NAEOffset                        : " << NAEOffset << endl;
    }

    // Compute dimensions of 4 sides in Mu2e coordinates
    double NAETopHalfX   = NAEHalfLengthXY + NAEHalfThickness;
    double NAETopHalfY   = NAEHalfThickness;
    double NAETopHalfZ   = NAEHalfLengthZ;
    double NAESideHalfX  = NAEHalfThickness;
    double NAESideHalfY  = NAEHalfLengthXY - NAEHalfThickness;
    double NAESideHalfZ  = NAEHalfLengthZ; 

    double NAETopDims[3] = {
      NAETopHalfX,
      NAETopHalfY,
      NAETopHalfZ
    };

    double NAESideDims[3] = {
      NAESideHalfX,
      NAESideHalfY,
      NAESideHalfZ
    };

    CLHEP::Hep3Vector NAETopOffset   (0.,  NAEHalfLengthXY, 0.);
    CLHEP::Hep3Vector NAEBottomOffset(0., -NAEHalfLengthXY, 0.);
    CLHEP::Hep3Vector NAELeftOffset  (     NAEHalfLengthXY, 0., 0.);
    CLHEP::Hep3Vector NAERightOffset (    -NAEHalfLengthXY, 0., 0.);


    VolumeInfo TopInfo    = nestBox("ExternalNeutronAbsorberTop",
                                    NAETopDims,
                                    NAMaterial,
                                    0,
                                    NAETopOffset + NAEOffset,
                                    hallInfo,
                                    0,
                                    NAVisible,
                                    G4Colour::Cyan(),
                                    NASolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );

    VolumeInfo BottomInfo = nestBox("ExternalNeutronAbsorberBottom",
                                    NAETopDims,
                                    NAMaterial,
                                    0,
                                    NAEBottomOffset + NAEOffset,
                                    hallInfo,
                                    0, 
                                    NAVisible,
                                    G4Colour::Cyan(), 
                                    NASolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );

    VolumeInfo LeftInfo   = nestBox("ExternalNeutronAbsorberLeft",
                                    NAESideDims,
                                    NAMaterial,
                                    0,
                                    NAELeftOffset + NAEOffset,
                                    hallInfo,
                                    0, 
                                    NAVisible,
                                    G4Colour::Cyan(),
                                    NASolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    ); 

    VolumeInfo RightInfo  = nestBox("ExternalNeutronAbsorberRight",
                                    NAESideDims,
                                    NAMaterial,
                                    0,
                                    NAERightOffset + NAEOffset,
                                    hallInfo,
                                    0, 
                                    NAVisible,
                                    G4Colour::Cyan(),
                                    NASolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    ); 

    // now constructing the internal neutron absorber
    // it is placed inside ToyDS2Vacuum & ToyDS3Vacuum like the protonabs1 & 2


    CLHEP::Hep3Vector NAI1OffsetInMu2e  = CLHEP::Hep3Vector(solenoidOffset,0.,NAIZ01);

    if ( verbosityLevel > 0) {
      cout << __func__ << " NAI1OffsetInMu2e                 : " << NAI1OffsetInMu2e << endl;
    }

    double NAIZ02 = NAIZ01 + NAIHalfLengthZ01 + NAIHalfLengthZ02;

    if ( verbosityLevel > 0) {
      cout << __func__ << " NAIZ01                           : " << NAIZ01  << endl;
      cout << __func__ << " NAIZ02                           : " << NAIZ02  << endl;
    }

    double NAI1FrontZ    = NAIZ01 - NAIHalfLengthZ01;
    double NAI12BoundaryZ = NAIZ01 + NAIHalfLengthZ01;
    double NAI2EndZ      = NAIZ02 + NAIHalfLengthZ02;

    if ( verbosityLevel > 0) {
      cout << __func__ << " NAI1FrontZ                       : " << NAI1FrontZ    << endl;
      cout << __func__ << " NAI12BoundaryZ                   : " << NAI12BoundaryZ << endl;
      cout << __func__ << " NAI2EndZ                         : " << NAI2EndZ      << endl;
    }

    // now local offsets in ToyDS2Vacuum & ToyDS3Vacuum

    VolumeInfo const & toyDS2VacuumInfo = _helper->locateVolInfo("ToyDS2Vacuum");
    VolumeInfo const & toyDS3VacuumInfo = _helper->locateVolInfo("ToyDS3Vacuum");

    // the Z offset of the InternalNeutronAbsorber2

    // we need to calculate where is the boundary of
    // InternalNeutronAbsorber1,2 wrt boudary between ToyDS2Vacuum &
    // ToyDS3Vacuum


    // G4VSolid* dsv1s = toyDS2VacuumInfo.solid
    // G4VPhysicalVolume* dsv1p = toyDS2VacuumInfo.physical;

    double ds2HalfLength = 
      dynamic_cast<G4Tubs* const>(toyDS2VacuumInfo.solid)->GetZHalfLength();

     double ds3HalfLength = 
      dynamic_cast<G4Tubs* const>(toyDS3VacuumInfo.solid->GetConstituentSolid(0))->GetZHalfLength();

   if ( verbosityLevel > 0) {
      cout << __func__ << " toyDS2VacuumInfo.solid->GetZHalfLength()    : " << 
        ds2HalfLength  << endl;

      cout << __func__ << " toyDS3VacuumInfo.solid->GetZHalfLength()    : " << 
        ds3HalfLength  << endl;

    }

    // what are the offsets of the DS2,3Vacuum in the Mu2e system?
    // the ds2/ds3 boundary in mu2e coordinates

    double ds23BoudaryZG = toyDS2VacuumInfo.centerInMu2e()(2) + 
      ds2HalfLength;

    double ds2FrontZG    = toyDS2VacuumInfo.centerInMu2e()(2) -
      ds2HalfLength;

    double ds3EndZG      = ds23BoudaryZG + 
      2.*ds3HalfLength;


    if ( verbosityLevel > 0) {
      cout << __func__ << " ds2FrontZG                       : " << ds2FrontZG    << endl;
      cout << __func__ << " ds23BoudaryZG                    : " << ds23BoudaryZG << endl;
      cout << __func__ << " ds3EndZG                         : " << ds3EndZG      << endl;
    }

    // certain combinations are illegal; we will assume that the
    // conical absorber is fully contained in ds2

    if ( ds2FrontZG>NAI1FrontZ ){
      throw cms::Exception("GEOM")
        << "Internal Neutron Absorber 1 front is outside of DS2Vaccum \n";
    }

    if ( ds23BoudaryZG<NAI12BoundaryZ ){
      throw cms::Exception("GEOM")
        << "Internal Neutron Absorber 1 end is outside of DS2Vaccum \n";
    }

    if ( ds3EndZG<NAI2EndZ ){
      throw cms::Exception("GEOM")
        << "Internal Neutron Absorber 2 end is outside of DS3Vaccum \n";
    }

    // tube parameters for InternalNeutronAbsorber2

    // we need to construct two InternalNeutronAbsorber2 one in ds2 and another one in ds3
    // we need to calculate their offsets


    // spliting the InternalNeutronAbsorber2

    
    double halfDeltaBoundary = (ds23BoudaryZG - NAI12BoundaryZ)*.5;
    double averageBoundary = (ds23BoudaryZG + NAI12BoundaryZ)*.5;

    CLHEP::Hep3Vector NAI23OffsetInMu2e = 
      CLHEP::Hep3Vector(solenoidOffset,0.,NAIZ02+halfDeltaBoundary);

    if ( verbosityLevel > 0) {
      cout << __func__ << " NAI23OffsetInMu2e                : " << NAI23OffsetInMu2e << endl;
    }

    CLHEP::Hep3Vector NAI22OffsetInMu2e = 
      CLHEP::Hep3Vector(solenoidOffset,0.,averageBoundary);

    if ( verbosityLevel > 0) {
      cout << __func__ << " NAI22OffsetInMu2e                : " << NAI22OffsetInMu2e << endl;
    }

    CLHEP::Hep3Vector NAI1Offset  =  NAI1OffsetInMu2e  - toyDS2VacuumInfo.centerInMu2e();
    CLHEP::Hep3Vector NAI22Offset =  NAI22OffsetInMu2e - toyDS2VacuumInfo.centerInMu2e();
    CLHEP::Hep3Vector NAI23Offset =  NAI23OffsetInMu2e - toyDS3VacuumInfo.centerInMu2e();

    if ( verbosityLevel > 0) {
      cout << __func__ << " toyDS2VacuumInfo.centerInMu2e()  : " << 
        toyDS2VacuumInfo.centerInMu2e() << endl;
      cout << __func__ << " toyDS3VacuumInfo.centerInMu2e()  : " << 
        toyDS3VacuumInfo.centerInMu2e() << endl;
      cout << __func__ << " NAI1Offset                       : " << NAI1Offset  << endl;
      cout << __func__ << " NAI22Offset                      : " << NAI22Offset << endl;
      cout << __func__ << " NAI23Offset                      : " << NAI23Offset << endl;
    }


    // constructing InternalNeutronAbsorber1 as an subtraction solid

    VolumeInfo INA1Info;
    INA1Info.name = "InternalNeutronAbsorber1";

    // full cylinder
    TubsParams NAI1Params( 0.,
                           NAIOuterRadius,
                           NAIHalfLengthZ01);

    G4Tubs* NAI1STubs = new G4Tubs(INA1Info.name + "Tubs", 
                                   NAI1Params.innerRadius,
                                   NAI1Params.outerRadius,
                                   NAI1Params.zHalfLength,
                                   NAI1Params.phi0, 
                                   NAI1Params.phiMax);

    // to be subtracted conical section
    G4Cons* NAI1SCons = new G4Cons(INA1Info.name + "Cons",
                                   0.,
                                   NAIInnerRadius0,
                                   0,
                                   NAIInnerRadius1,
                                   NAIHalfLengthZ01,
                                   0.,
                                   CLHEP::twopi);

    INA1Info.solid = new G4SubtractionSolid(INA1Info.name, NAI1STubs, NAI1SCons);

    finishNesting(INA1Info,
                  NAMaterial,
                  0,
                  NAI1Offset,
                  toyDS2VacuumInfo.logical,
                  0,
                  NAVisible,
                  G4Colour::Cyan(),
                  NASolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    // the "split" cylindrical internal part

    TubsParams NAI22Params( NAIInnerRadius2,
                            NAIOuterRadius,
                            halfDeltaBoundary);

    TubsParams NAI23Params( NAIInnerRadius2,
                            NAIOuterRadius,
                            NAIHalfLengthZ02-halfDeltaBoundary);

 
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
