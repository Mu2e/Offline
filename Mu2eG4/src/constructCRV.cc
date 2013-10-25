//
// Free function to create CRV aka Scintillator Shield in CosmicRayShield
//
// $Id: constructCRV.cc,v 1.19 2013/10/25 05:06:33 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/10/25 05:06:33 $
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
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"
#include "Mu2eG4/inc/nestBox.hh"

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

namespace mu2e 
{
  void constructCRV( VolumeInfo const & parent, SimpleConfig const &  _config)
  {
    GeomHandle<CosmicRayShield> CosmicRayShieldGeomHandle;

    G4Helper    & _helper = *(art::ServiceHandle<G4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    bool scintillatorShieldVisible = _config.getBool("crs.vetoVisible",true);
    bool scintillatorShieldSolid   = _config.getBool("crs.vetoSolid",false);

    int verbosityLevel = _config.getInt("crs.verbosityLevel",0);

    bool const forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false);

    std::string hallAirMaterialName = _config.getString("hall.insideMaterialName");

    std::map<std::string,CRSScintillatorShield> const &shields = CosmicRayShieldGeomHandle->getCRSScintillatorShields();

    CLHEP::Hep3Vector parentCenterInMu2e = parent.centerInMu2e();

    std::map<std::string,CRSScintillatorShield>::const_iterator ishield;
    for(ishield=shields.begin(); ishield!=shields.end(); ++ishield) 
    {
      CRSScintillatorShield const & shield = ishield->second;
      std::string shieldName = ishield->first;

      if(verbosityLevel > 0) cout << __func__ << " constructing            : " << shieldName << endl;

      // all materials and dimensions are the same within a particular shield
      CRSScintillatorBarDetail const & barDetail = shield.getCRSScintillatorBarDetail();

      G4Material* scintillatorBarMaterial = findMaterialOrThrow(barDetail.getMaterialName());
      std::vector<double> const &  scintillatorBarHalfLengths = barDetail.getHalfLengths();

      std::string scintillatorBarName = shieldName;

      G4VSolid* scintillatorBarSolid = new G4Box(scintillatorBarName,
                                                 scintillatorBarHalfLengths[0],
                                                 scintillatorBarHalfLengths[1],
                                                 scintillatorBarHalfLengths[2]);

      G4LogicalVolume* scintillatorBarLogical = new G4LogicalVolume(scintillatorBarSolid,
                                                                    scintillatorBarMaterial,
                                                                    scintillatorBarName);

      // visibility attributes
      if (!scintillatorShieldVisible) 
      {
        scintillatorBarLogical->SetVisAttributes(G4VisAttributes::Invisible);
      }
      else 
      {
        G4Colour  orange  (.75, .55, .0);
        G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, orange));
        visAtt->SetForceSolid(scintillatorShieldSolid);
        visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
        scintillatorBarLogical->SetVisAttributes(visAtt);
      }

      // Make each scintillatorBar a sensitive detector.
      G4VSensitiveDetector *sd = 
                     G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CRSScintillatorBar());
      if(sd) scintillatorBarLogical->SetSensitiveDetector(sd);

      if (verbosityLevel > 1) 
      {
        G4SDManager::GetSDMpointer()->SetVerboseLevel(verbosityLevel-1);
      }

      const std::vector<CRSScintillatorModule> &modules = shield.getCRSScintillatorModules();
      std::vector<CRSScintillatorModule>::const_iterator imodule;
      for(imodule=modules.begin(); imodule!=modules.end(); ++imodule) 
      {
        const std::vector<CRSScintillatorLayer> &layers = imodule->getLayers();
        std::vector<CRSScintillatorLayer>::const_iterator ilayer;
        for(ilayer=layers.begin(); ilayer!=layers.end(); ++ilayer) 
        {

          //construct a G4Box around all bars of a layer to speed up the surface checks
          std::vector<double> layerHalflengths(3);
          CLHEP::Hep3Vector   layerCenterInMu2e;
          ilayer->getDimensions(layerHalflengths, layerCenterInMu2e);
          CLHEP::Hep3Vector layerAirOffset = layerCenterInMu2e - parentCenterInMu2e;

          G4Material* hallAirMaterial = findMaterialOrThrow(hallAirMaterialName);
          VolumeInfo layerInfo = nestBox(ilayer->name("CRSScintillatorLayer_"),
                                         layerHalflengths,
                                         hallAirMaterial,
                                         NULL/*rotation*/,
                                         layerAirOffset,
                                         parent.logical,
                                         0/*copyNumber*/,
                                         false/*visible*/,
                                         G4Colour::Yellow(),
                                         false/*solid*/,
                                         forceAuxEdgeVisible,
                                         true/*placePV*/,
                                         doSurfaceCheck);

          const std::vector<const CRSScintillatorBar*> &bars = ilayer->getBars();
          std::vector<const CRSScintillatorBar*>::const_iterator ibar;
          for(ibar=bars.begin(); ibar!=bars.end(); ++ibar) 
          {
            const CRSScintillatorBar &bar = **ibar; 

            //bar.getPosition() returns the bar position in Mu2e coordinates
            CLHEP::Hep3Vector barLayerOffset = bar.getPosition() - layerCenterInMu2e; 

            if ( verbosityLevel > 3 ) 
            {
              cout << __func__ << " barPosition          : " <<  bar.getPosition() << endl;
              cout << __func__ << " barShieldOffset      : " <<  barLayerOffset << endl;
              cout << __func__ << " shieldAirOffset      : " <<  layerAirOffset << endl;
            }

            G4VPhysicalVolume* pv = new G4PVPlacement( NULL,
                                                       barLayerOffset,
                                                       scintillatorBarLogical,
                                                       bar.name("CRSScintillatorBar_"),
                                                       layerInfo.logical,
                                                       false,
                                                       bar.index().asInt(),
                                                       false);

            if ( doSurfaceCheck) 
            {
              checkForOverlaps( pv, _config, verbosityLevel>0);
            }
          } //ibar
        } //ilayer
      } //imodule
    } //ishield
  } //construct CRV

} //namespace mu2e
