//
// Free function to create Neutron Absorbers in G4
//
// $Id: constructExternalNeutronAbsorber.cc,v 1.2 2013/06/12 19:52:03 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/06/12 19:52:03 $
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
#include "ExternalNeutronAbsorberGeom/inc/ExternalNeutronAbsorber.hh"
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

    // Flags
    int  const verbosityLevel      = _config.getInt("extneutronabs.verbosityLevel",0);
    bool const NAVisible           = _config.getBool("extneutronabs.visible",false);
    bool const NASolid             = _config.getBool("extneutronabs.solid");
    bool const forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    // Fetch ext. neutron absorber geometry
    GeomHandle<ExternalNeutronAbsorber> nae;
    G4Material* NAMaterial = findMaterialOrThrow( nae->material() );

    // Access to the G4HelperService - get hall volume
    G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));
    VolumeInfo const & hallInfo = _helper->locateVolInfo("HallAir");

    if ( verbosityLevel > 0) {
      cout << __func__ << " NAEOffsetInMu2e                  : " << nae->position().z() << endl;
    }
    
    // now local offset in hallAir
    CLHEP::Hep3Vector NAEOffset =  nae->position() - hallInfo.centerInMu2e();
    
    if ( verbosityLevel > 0) {
      cout << __func__ << " hallInfo.centerInMu2e()          : " << hallInfo.centerInMu2e() << endl;
      cout << __func__ << " NAEOffset                        : " << NAEOffset << endl;
    }
    
    // Compute dimensions of 4 sides in Mu2e coordinates
    double NAETopHalfX   = nae->halfLengthXY() + nae->halfThickness();
    double NAETopHalfY   = nae->halfThickness();
    double NAETopHalfZ   = nae->halfLengthZ();
    double NAESideHalfX  = nae->halfThickness();
    double NAESideHalfY  = nae->halfLengthXY() - nae->halfThickness();
    double NAESideHalfZ  = nae->halfLengthZ();
    
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
    
    CLHEP::Hep3Vector NAETopOffset   (0.,  nae->halfLengthXY(), 0.);
    CLHEP::Hep3Vector NAEBottomOffset(0., -nae->halfLengthXY(), 0.);
    CLHEP::Hep3Vector NAELeftOffset  (     nae->halfLengthXY(), 0., 0.);
    CLHEP::Hep3Vector NAERightOffset (    -nae->halfLengthXY(), 0., 0.);
    
    
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


