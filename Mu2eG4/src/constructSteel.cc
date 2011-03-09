//
// Free function to create Hall Steel
//
// $Id: constructSteel.cc,v 1.3 2011/03/09 19:51:42 genser Exp $
// $Author: genser $
// $Date: 2011/03/09 19:51:42 $
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
#include "CosmicRayShieldGeom/inc/CRSSteelShield.hh"
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

    std::map<std::string,CRSSteelShield> const & shields = 
      CosmicRayShieldGeomHandle->getCRSSteelShields();

    for (std::map<std::string,CRSSteelShield>::const_iterator ishield=shields.begin();
         ishield!=shields.end(); ++ishield) {

      CRSSteelShield const & shield = ishield->second;
      std::string shieldName = ishield->first;

      if (shield.getHoleRadius() == 0.) {
        nestBox(shield.name(),
                shield.getHalfLengths(),
                CRSSteelShieldMaterial,
                shield.getRotation(),
                shield.getLocalOffset(),
                parent,
                0,
                CRSSteelShieldVisible,
                G4Colour::Green(),
                CRSSteelShieldSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );      
      } else {

        // constructing "hollow" upstream steel

        G4Box* CRSSteelShieldBox = new G4Box(shield.name()+"Box",
                                  shield.getHalfLengths()[0], 
                                  shield.getHalfLengths()[1], 
                                  shield.getHalfLengths()[2]);

        //Hole in upstream shield for TS is centered in the shield

        TubsParams UpstreamHoleDims(0.,
                                 shield.getHoleRadius(),
                                 shield.getHalfLengths()[2]);

        G4Tubs* CRSSteelShieldUpstreamHoleTubs = new G4Tubs(shield.name()+"HoleTubs", 
                                                    UpstreamHoleDims.innerRadius,
                                                    UpstreamHoleDims.outerRadius,
                                                    UpstreamHoleDims.zHalfLength,
                                                    UpstreamHoleDims.phi0, 
                                                    UpstreamHoleDims.phiMax);

        VolumeInfo CRSSteelShieldUpstreamInfo;

        CRSSteelShieldUpstreamInfo.name = shield.name();

        CRSSteelShieldUpstreamInfo.solid = 
          new G4SubtractionSolid(CRSSteelShieldUpstreamInfo.name, 
                                 CRSSteelShieldBox, CRSSteelShieldUpstreamHoleTubs);

        finishNesting(CRSSteelShieldUpstreamInfo,
                      CRSSteelShieldMaterial,
                      shield.getRotation(),
                      shield.getLocalOffset(),
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

  }

}
