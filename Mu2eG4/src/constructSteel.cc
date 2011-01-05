//
// Free function to create Hall Steel
//
// $Id: constructSteel.cc,v 1.1 2011/01/05 21:04:47 genser Exp $
// $Author: genser $
// $Date: 2011/01/05 21:04:47 $
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

    // Extract information from the config file
    double HallSteelHalfThick   = _config->getDouble("fluxcrv.HallSteelHalfThick");
    double HallSteelHalfLenXY = _config->getDouble("fluxcrv.HallSteelHalfLengthXY");
    double HallSteelHalfLenZ = _config->getDouble("fluxcrv.HallSteelHalfLengthZ");
    G4Material* HallSteelShieldMaterial = materialFinder.get("fluxcrv.HallSteelMaterialName");

    GeomHandle<Beamline> beamg;

    // Compute dimensions of 5 sides in Mu2e coordinates
    double HallSteelTopHalfX   = HallSteelHalfLenXY + HallSteelHalfThick;
    double HallSteelTopHalfY   = HallSteelHalfThick;
    double HallSteelTopHalfZ   = HallSteelHalfLenZ;
    double HallSteelSideHalfX  = HallSteelHalfThick;
    double HallSteelSideHalfY  = HallSteelHalfLenXY - HallSteelHalfThick;
    double HallSteelSideHalfZ  = HallSteelHalfLenZ; 
    double HallSteelFrontHalfX = HallSteelHalfLenXY - HallSteelHalfThick;
    double HallSteelFrontHalfY = HallSteelHalfLenXY - HallSteelHalfThick;
    double HallSteelFrontHalfZ = HallSteelHalfThick;

    double HallSteelTopDims[3] ={
      HallSteelTopHalfX,
      HallSteelTopHalfY,
      HallSteelTopHalfZ
    };

    double HallSteelSideDims[3] ={
      HallSteelSideHalfX,
      HallSteelSideHalfY,
      HallSteelSideHalfZ
    };

    double HallSteelFrontDims[3] ={
      HallSteelFrontHalfX,
      HallSteelFrontHalfY,           
      HallSteelFrontHalfZ                        
    };
    TubsParams FrontHoleDims(0.,
                             _config->getDouble("fluxcrv.HallSteelHoleRadius"),
                             HallSteelHalfThick
                             );

    // Get positions of each side. Assuming view from target foils 
    double dsCoilZ0          = _config->getDouble("toyDS.z0");

    double solenoidOffset = beamg->solenoidOffset();

    //G4ThreeVector detSolCoilPosition(-solenoidOffset, 0., -dsCoilZ0);
    G4ThreeVector detSolCoilPosition(+solenoidOffset, 0., -dsCoilZ0);

    std::vector<double> HallSteelOffsetSTDV;
    _config->getVectorDouble("fluxcrv.HallSteelOffset", HallSteelOffsetSTDV, 3);
    G4ThreeVector HallSteelOffset(HallSteelOffsetSTDV[0],
                                  HallSteelOffsetSTDV[1],
                                  HallSteelOffsetSTDV[2]);

    // Imagine a box that exactly contains the flux return steel.
    // This is the center of that box, in the coordinate system of the mother volume(the hall air).
    CLHEP::Hep3Vector boxCenter = -(parent.centerInWorld - VolumeInfo::getMu2eOriginInWorld() 
                                    + detSolCoilPosition + HallSteelOffset);

    G4ThreeVector TopShield   (0.,      HallSteelSideHalfY + HallSteelHalfThick, 0.);
    G4ThreeVector BottomShield(0.,    -(HallSteelSideHalfY + HallSteelHalfThick), 0.);
    G4ThreeVector LeftShield  (         HallSteelSideHalfY + HallSteelHalfThick,0., 0.);
    G4ThreeVector RightShield (       -(HallSteelSideHalfY + HallSteelHalfThick),0., 0.);
    G4ThreeVector BackShield  (0., 0.,  HallSteelSideHalfZ - HallSteelHalfThick);
    G4ThreeVector FrontShield (0., 0.,-(HallSteelSideHalfZ - HallSteelHalfThick));

    //Hole in front shield for TS is centered in the shield
    G4ThreeVector FrontHole(0.,0.,0.);

    bool hallSteelVisible = _config->getBool("fluxcrv.visible",true);
    bool hallSteelSolid   = _config->getBool("fluxcrv.solid",false);

    bool const forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    // Place Boxes

    VolumeInfo TopInfo = nestBox ("HallSteelTopShield",
                                  HallSteelTopDims,
                                  HallSteelShieldMaterial,
                                  0,
                                  TopShield + boxCenter,
                                  parent,
                                  0,
                                  hallSteelVisible,
                                  G4Colour::Green(),
                                  hallSteelSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck
                                  );

    VolumeInfo BottomInfo = nestBox ("HallSteelBottomShield",
                                     HallSteelTopDims, 
                                     HallSteelShieldMaterial,
                                     0,
                                     BottomShield + boxCenter,
                                     parent,
                                     0, 
                                     hallSteelVisible,
                                     G4Colour::Green(), 
                                     hallSteelSolid,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     doSurfaceCheck
                                     );

    VolumeInfo LeftInfo = nestBox ("HallSteelLeftShield",
                                   HallSteelSideDims,
                                   HallSteelShieldMaterial,
                                   0, 
                                   LeftShield + boxCenter,
                                   parent,
                                   0, 
                                   hallSteelVisible,
                                   G4Colour::Green(),
                                   hallSteelSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ); 

    VolumeInfo RightInfo = nestBox ("HallSteelRightShield",
                                    HallSteelSideDims,
                                    HallSteelShieldMaterial,
                                    0, 
                                    RightShield + boxCenter,
                                    parent,
                                    0, 
                                    hallSteelVisible,
                                    G4Colour::Green(),
                                    hallSteelSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    ); 

    VolumeInfo BackInfo = nestBox ("HallSteelBackShield",
                                   HallSteelFrontDims,
                                   HallSteelShieldMaterial,
                                   0, 
                                   BackShield + boxCenter,
                                   parent,
                                   0, 
                                   hallSteelVisible,
                                   G4Colour::Green(),
                                   hallSteelSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ); 

    // constructing "hollow" front shield

    G4Box* CRVBox = new G4Box("CRVFrontShieldBox", 
                              HallSteelFrontDims[0], 
                              HallSteelFrontDims[1], 
                              HallSteelFrontDims[2]);

    G4Tubs* HallSteelFrontHoleTubs = new G4Tubs("HallSteelFrontHoleTubs", 
                                                FrontHoleDims.innerRadius,
                                                FrontHoleDims.outerRadius,
                                                FrontHoleDims.zHalfLength,
                                                FrontHoleDims.phi0, 
                                                FrontHoleDims.phiMax);

    VolumeInfo FrontShieldInfo;

    FrontShieldInfo.name = "CRVFrontShield";

    FrontShieldInfo.solid = 
      new G4SubtractionSolid(FrontShieldInfo.name, CRVBox, HallSteelFrontHoleTubs);

    finishNesting(FrontShieldInfo,
                  HallSteelShieldMaterial,
                  0,
                  FrontShield  + boxCenter, 
                  parent.logical,
                  0,
                  hallSteelVisible,
                  G4Colour::Green(),
                  hallSteelSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

  }


}
