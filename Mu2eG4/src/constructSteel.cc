//
// Free function to create Hall Steel
//
// $Id: constructSteel.cc,v 1.2 2011/01/25 16:47:55 genser Exp $
// $Author: genser $
// $Date: 2011/01/25 16:47:55 $
//
// Original author KLG based on Mu2eWorld constructSteel
//
// Notes:

// Mu2e includes.

#include "Mu2eG4/inc/constructSteel.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShieldSteelShield.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "TrackerGeom/inc/TubsParams.hh"


// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"


using namespace std;

namespace mu2e {

  void constructSteel( const VolumeInfo& parent, 
                       SimpleConfig const * const _config
                       ){

    MaterialFinder materialFinder(*_config);

    G4Material* CRSSteelShieldMaterial = materialFinder.get("fluxcrv.HallSteelMaterialName");

    bool CRSSteelShieldVisible = _config->getBool("fluxcrv.visible",true);
    bool CRSSteelShieldSolid   = _config->getBool("fluxcrv.solid",false);

    bool const forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    // get the CRS parameters from the geometry service and place the steel boxes

    GeomHandle<CosmicRayShield> CosmicRayShieldGeomHandle;

    CosmicRayShieldSteelShield const & CRSSteelTopShield = 
      CosmicRayShieldGeomHandle->getCosmicRayShieldSteelShield("CRSSteelTopShield");

    VolumeInfo TopInfo   = nestBox(CRSSteelTopShield.name(),
                                   CRSSteelTopShield.getHalfLengths(),
                                   CRSSteelShieldMaterial,
                                   CRSSteelTopShield.getRotation(),
                                   CRSSteelTopShield.getLocalOffset(),
                                   parent,
                                   0,
                                   CRSSteelShieldVisible,
                                   G4Colour::Green(),
                                   CRSSteelShieldSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   );

    CosmicRayShieldSteelShield const & CRSSteelBottomShield = 
      CosmicRayShieldGeomHandle->getCosmicRayShieldSteelShield("CRSSteelBottomShield");

    VolumeInfo BottomInfo = nestBox(CRSSteelBottomShield.name(),
                                    CRSSteelBottomShield.getHalfLengths(),
                                    CRSSteelShieldMaterial,
                                    CRSSteelBottomShield.getRotation(),
                                    CRSSteelBottomShield.getLocalOffset(),
                                    parent,
                                    0, 
                                    CRSSteelShieldVisible,
                                    G4Colour::Green(), 
                                    CRSSteelShieldSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );

    CosmicRayShieldSteelShield const & CRSSteelLeftShield = 
      CosmicRayShieldGeomHandle->getCosmicRayShieldSteelShield("CRSSteelLeftShield");

    VolumeInfo LeftInfo   = nestBox(CRSSteelLeftShield.name(),
                                    CRSSteelLeftShield.getHalfLengths(),
                                    CRSSteelShieldMaterial,
                                    CRSSteelLeftShield.getRotation(),
                                    CRSSteelLeftShield.getLocalOffset(),
                                    parent,
                                    0, 
                                    CRSSteelShieldVisible,
                                    G4Colour::Green(),
                                    CRSSteelShieldSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    ); 

    CosmicRayShieldSteelShield const & CRSSteelRightShield = 
      CosmicRayShieldGeomHandle->getCosmicRayShieldSteelShield("CRSSteelRightShield");

    VolumeInfo RightInfo = nestBox(CRSSteelRightShield.name(),
                                   CRSSteelRightShield.getHalfLengths(),
                                   CRSSteelShieldMaterial,
                                   CRSSteelRightShield.getRotation(),
                                   CRSSteelRightShield.getLocalOffset(),
                                   parent,
                                   0, 
                                   CRSSteelShieldVisible,
                                   G4Colour::Green(),
                                   CRSSteelShieldSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ); 

    CosmicRayShieldSteelShield const & CRSSteelBackShield = 
      CosmicRayShieldGeomHandle->getCosmicRayShieldSteelShield("CRSSteelBackShield");

    VolumeInfo BackInfo = nestBox(CRSSteelBackShield.name(),
                                  CRSSteelBackShield.getHalfLengths(),
                                  CRSSteelShieldMaterial,
                                  CRSSteelBackShield.getRotation(),
                                  CRSSteelBackShield.getLocalOffset(),
                                  parent,
                                  0, 
                                  CRSSteelShieldVisible,
                                  G4Colour::Green(),
                                  CRSSteelShieldSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck
                                  ); 

    // constructing "hollow" front shield

    CosmicRayShieldSteelShield const & CRSSteelFrontShield = 
      CosmicRayShieldGeomHandle->getCosmicRayShieldSteelShield("CRSSteelFrontShield");


    G4Box* CRVBox = new G4Box(CRSSteelFrontShield.name()+"Box",
                              CRSSteelFrontShield.getHalfLengths()[0], 
                              CRSSteelFrontShield.getHalfLengths()[1], 
                              CRSSteelFrontShield.getHalfLengths()[2]);

    //Hole in front shield for TS is centered in the shield

    TubsParams FrontHoleDims(0.,
                             CRSSteelFrontShield.getHoleRadius(),
                             CRSSteelFrontShield.getHalfLengths()[2]
                             );

    G4Tubs* HallSteelFrontHoleTubs = new G4Tubs(CRSSteelFrontShield.name()+"HoleTubs", 
                                                FrontHoleDims.innerRadius,
                                                FrontHoleDims.outerRadius,
                                                FrontHoleDims.zHalfLength,
                                                FrontHoleDims.phi0, 
                                                FrontHoleDims.phiMax);

    VolumeInfo FrontShieldInfo;

    FrontShieldInfo.name = CRSSteelFrontShield.name();

    FrontShieldInfo.solid = 
      new G4SubtractionSolid(FrontShieldInfo.name, CRVBox, HallSteelFrontHoleTubs);

    finishNesting(FrontShieldInfo,
                  CRSSteelShieldMaterial,
                  CRSSteelFrontShield.getRotation(),
                  CRSSteelFrontShield.getLocalOffset(),
                  parent.logical,
                  0,
                  CRSSteelShieldVisible,
                  G4Colour::Green(),
                  CRSSteelShieldSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

  }


}
