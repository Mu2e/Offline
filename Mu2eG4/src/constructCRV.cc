//
// Free function to create CRV aka Scintillator Shield in CosmicRayShield
//
// $Id: constructCRV.cc,v 1.5 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:18 $
//
// Original author KLG
//
// Notes:


// clhep includes
#include "CLHEP/Vector/ThreeVector.h"

// Mu2e includes.

#include "Mu2eG4/inc/constructCRV.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"

#include "Mu2eG4/inc/CRSScintillatorBarSD.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"

// G4 includes

#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Box.hh"

#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4VisAttributes.hh"

#include "G4RotationMatrix.hh"

#include "G4SDManager.hh"

using namespace std;

namespace mu2e {

  void constructCRV( const VolumeInfo& parent,
                     SimpleConfig const * const _config
                     ){

    G4Helper    & _helper = *(art::ServiceHandle<G4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    bool scintillatorShieldVisible = _config->get<bool>("crs.vetoVisible",true);
    bool scintillatorShieldSolid   = _config->get<bool>("crs.vetoSolid",false);

    int verbosityLevel = _config->get<int>("crs.verbosityLevel",0);

    bool const forceAuxEdgeVisible = _config->get<bool>("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config->get<bool>("g4.doSurfaceCheck",false);

    // get the CRS parameters from the geometry service and place the veto elements

    GeomHandle<CosmicRayShield> CosmicRayShieldGeomHandle;

    std::map<std::string,CRSScintillatorShield> const & shields =
      CosmicRayShieldGeomHandle->getCRSScintillatorShields();

    CLHEP::Hep3Vector perentCenterInMu2e = parent.centerInWorld - VolumeInfo::getMu2eOriginInWorld();

    // all materials and dimensions are the same

    CRSScintillatorBarDetail const & barDetail =
      CosmicRayShieldGeomHandle->getCRSScintillatorBarDetail();

    G4Material* scintillatorBarMaterial =
      findMaterialOrThrow(barDetail.getMaterialName(0));

    std::vector<double> const &  scintillatorBarHalfLengths = barDetail.getHalfLengths();

    // each solid is the same,
    // is each logical volume the same?
    // but each physical volume has different name...
    // this seems not quite compatible with VolumeInfo and their registry, so we will not use it

    //    VolumeInfo scintillatorBarInfo;

    std::string scintillatorBarName = "CRSScintillatorBar";

    G4VSolid* scintillatorBarSolid =
      new G4Box( scintillatorBarName,
                 scintillatorBarHalfLengths[0],
                 scintillatorBarHalfLengths[1],
                 scintillatorBarHalfLengths[2]);

    G4LogicalVolume* scintillatorBarLogical =
      new G4LogicalVolume( scintillatorBarSolid,
                           scintillatorBarMaterial,
                           scintillatorBarName);

    // visibility attributes

    if (!scintillatorShieldVisible) {

      scintillatorBarLogical->SetVisAttributes(G4VisAttributes::Invisible);

    } else {
      G4Colour  orange  (.75, .55, .0);
      G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, orange));
      visAtt->SetForceSolid(scintillatorShieldSolid);
      visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
      scintillatorBarLogical->SetVisAttributes(visAtt);

    }

    // Make each scintillatorBar a sensitive detector.

    scintillatorBarLogical->
      SetSensitiveDetector(G4SDManager::GetSDMpointer()->
                           FindSensitiveDetector(SensitiveDetectorName::CRSScintillatorBar()) );

    if (verbosityLevel > 1) {
      G4SDManager::GetSDMpointer()->SetVerboseLevel(verbosityLevel-1);
    }

    for (std::map<std::string,CRSScintillatorShield>::const_iterator ishield=shields.begin();
         ishield!=shields.end(); ++ishield) {

      CRSScintillatorShield const & shield = ishield->second;
      std::string shieldName = ishield->first;

      //      if (shieldName!="CRSScintillatorDShield") continue;

      verbosityLevel > 0 &&
        cout << __func__ << " constructing            : " << shieldName << endl;

      // we shall place individual scintillators skipping modules and layers for now
      // the scintillators do have some common parameters though

      CLHEP::HepRotationX RX(shield.getGlobalRotationAngles()[0]);
      CLHEP::HepRotationY RY(shield.getGlobalRotationAngles()[1]);
      CLHEP::HepRotationZ RZ(shield.getGlobalRotationAngles()[2]);

      G4RotationMatrix* shieldRotation = reg.add(G4RotationMatrix(RX*RY*RZ));

      if ( verbosityLevel > 0 ) {
        cout << __func__ << " shieldLocalOffset       : " << shield.getLocalOffset() << endl;
        cout << __func__ << " shieldGlobalOffset      : " << shield.getGlobalOffset() << endl;
        cout << __func__ << " shieldAirOffset         : " << shield.getGlobalOffset() - perentCenterInMu2e << endl;
        cout << __func__ << " getGlobalRotationAngles : " <<
          shield.getGlobalRotationAngles()[0] << ", " <<
          shield.getGlobalRotationAngles()[1] << ", " <<
          shield.getGlobalRotationAngles()[2] << ", " << endl;
        cout << __func__ << " shieldRotation and *    : " << shieldRotation << *shieldRotation << endl;
      }

      // now loop over all bars in the given shield

      int nModules = shield.nModules();

      verbosityLevel > 0 &&
        cout << __func__ << " nModules                : " << nModules << endl;

      for (int im = 0; im < nModules; ++im) {

        verbosityLevel > 0 &&
          cout << __func__ << " working on module       : " << im << endl;

        CRSScintillatorModule const & module = shield.getModule(im);

        if ( verbosityLevel > 1 ) {
          cout << __func__ << " moduleLocalOffset       : " << module.getLocalOffset() << endl;
          cout << __func__ << " moduleGlobalOffset      : " << module.getGlobalOffset() << endl;
          cout << __func__ << " moduleAirOffset         : " << module.getGlobalOffset() - perentCenterInMu2e << endl;
        }

        int nLayers = module.nLayers();

        verbosityLevel > 0 &&
          cout << __func__ << " nLayers                 : " << nLayers << endl;

        for (int il = 0; il < nLayers; ++il) {

          verbosityLevel > 1 &&
            cout << __func__ << " working on layer        : " << il << endl;

          CRSScintillatorLayer const & layer = module.getLayer(il);

          if ( verbosityLevel > 2 ) {
            cout << __func__ << " layerLocalOffset        : " << layer.getLocalOffset() << endl;
            cout << __func__ << " layerGlobalOffset       : " << layer.getGlobalOffset() << endl;
            cout << __func__ << " layerAirOffset          : " << layer.getGlobalOffset() - perentCenterInMu2e << endl;
          }

          int nBars = layer.nBars();

          verbosityLevel > 0 &&
            cout << __func__ << " nBars                   : " << nBars << endl;

          for (int ib = 0; ib < nBars; ++ib) {

            verbosityLevel > 2 &&
              cout << __func__ << " working on bar          : " << ib << endl;

            CRSScintillatorBar const & bar = layer.getBar(ib);

            // placing the bar;

            // we need the offsets "local" to the air
            // so we need the globaloffsets for each bar and the air (parent volume)

            CLHEP::Hep3Vector barAirOffset = bar.getGlobalOffset() - perentCenterInMu2e;

            // the rotation is the same for each bar in a given shield

            if ( verbosityLevel > 3 ) {
              cout << __func__ << " barLocalOffset       : " <<  bar.getLocalOffset() << endl;
              cout << __func__ << " barGlobalOffset      : " <<  bar.getGlobalOffset() << endl;
              cout << __func__ << " barAirOffset         : " <<  barAirOffset << endl;
            }

            new G4PVPlacement( shieldRotation,
                               barAirOffset,
                               scintillatorBarLogical,
                               bar.name(scintillatorBarName+"_"),
                               parent.logical,
                               0,
                               bar.Index().asInt(),
                               doSurfaceCheck);

          }

        }

      }

    }

  }

}
