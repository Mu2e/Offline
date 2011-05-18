//
// Free function to create Hall Steel
//
// $Id: constructSteel.cc,v 1.6 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:18 $
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

    int const verbosityLevel = _config->get<int>("crs.verbosityLevel",0);

    MaterialFinder materialFinder(*_config);

    G4Material* CRSSteelShieldMaterial = materialFinder.get("fluxcrv.HallSteelMaterialName");

    bool CRSSteelShieldVisible = _config->get<bool>("fluxcrv.visible",true);
    bool CRSSteelShieldSolid   = _config->get<bool>("fluxcrv.solid",false);

    bool const forceAuxEdgeVisible = _config->get<bool>("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config->get<bool>("g4.doSurfaceCheck",false);
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

        // constructing "hollow" upstream/downstream steel

        G4Box* CRSSteelShieldBox = new G4Box(shield.name()+"Box",
                                  shield.getHalfLengths()[0],
                                  shield.getHalfLengths()[1],
                                  shield.getHalfLengths()[2]);

        //Hole in  shield for TS is centered in the shield

        TubsParams HoleDims(0.,
                            shield.getHoleRadius(),
                            shield.getHalfLengths()[2]+1.0);//  added for better rendering

        G4Tubs* CRSSteelShieldHoleTubs = new G4Tubs(shield.name()+"HoleTubs",
                                                    HoleDims.innerRadius,
                                                    HoleDims.outerRadius,
                                                    HoleDims.zHalfLength,
                                                    HoleDims.phi0,
                                                    HoleDims.phiMax);

        VolumeInfo CRSSteelShieldInfo;

        CRSSteelShieldInfo.name = shield.name();

        CRSSteelShieldInfo.solid =
          new G4SubtractionSolid(CRSSteelShieldInfo.name,
                                 CRSSteelShieldBox, CRSSteelShieldHoleTubs);

        finishNesting(CRSSteelShieldInfo,
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

        if ( verbosityLevel > 0) {
          double zhl  =  shield.getHalfLengths()[2];
          cout << __func__ << " " << shield.name() << " zhl    : " << zhl << endl;
         double CRSSteelShieldOffsetInMu2eZ = shield.getGlobalOffset()[CLHEP::Hep3Vector::Z];
         cout << __func__ << " " << shield.name() << " Z extent in Mu2e    : " <<
            CRSSteelShieldOffsetInMu2eZ - zhl << ", " << CRSSteelShieldOffsetInMu2eZ + zhl << endl;
        }

      }

    }

  }

}
