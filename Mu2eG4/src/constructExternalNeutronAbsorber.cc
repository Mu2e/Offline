//
// Free function to create Neutron Absorbers in G4
//
// $Id: constructExternalNeutronAbsorber.cc,v 1.1 2013/05/31 15:58:59 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/05/31 15:58:59 $
//
// Original author KLG
//
// Notes:
// Construct the Neutron Absorbers in G4


// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Mu2e includes.

#include "Mu2eG4/inc/constructExternalNeutronAbsorber.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
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
#include "G4Polycone.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4VPhysicalVolume.hh"

using namespace std;

namespace mu2e {

  void constructExternalNeutronAbsorber(SimpleConfig const & _config){

    int const verbosityLevel = _config.getInt("neutronabsorber.verbosityLevel",0);

    // the absorber is split into two major pieces internal & external (wrt to the coil)
    // the external part is placed in the hall air (and does not need to be split)

    // Extract information from the config file.

    string NAMaterialName      = _config.getString("neutronabsorber.materialName");
    double NAEHalfLengthZ      = _config.getDouble("neutronabsorber.externalHalfLengthZ");
    double NAEHalfLengthXY     = _config.getDouble("neutronabsorber.externalHalfLengthXY");
    double NAEHalfThickness    = _config.getDouble("neutronabsorber.externalHalfThickness");
    double NAEZ0               = _config.getDouble("neutronabsorber.externalZ0");

    bool NAVisible             = _config.getBool("neutronabsorber.visible",false);
    bool NASolid               = _config.getBool("neutronabsorber.solid");

    bool const forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    G4Material* NAMaterial = findMaterialOrThrow(NAMaterialName);

    // constructing External Neutron Absorber, placing it in HallAir (the name is hardcoded here...)

    // Access to the G4HelperService.
    G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));

    VolumeInfo const & hallInfo = _helper->locateVolInfo("HallAir");

    // NAEZ0 is the offset in mu2e
    // the Z offset in World is NAEZ0+mu2eOrigin

    // Get DS geom object
    GeomHandle<DetectorSolenoid> ds;

    if ( verbosityLevel > 0) {
      cout << __func__ << " solenoidOffset                   : " << ds->position().x()  << endl;
    }

    CLHEP::Hep3Vector NAEOffsetInMu2e  = CLHEP::Hep3Vector(ds->position().x(),0.,NAEZ0);
    
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


  }

} // end of constructExternalNeutronAbsorber;


