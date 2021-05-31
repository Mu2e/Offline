//
// Free function to create CRV aka Scintillator Shield in CosmicRayShield
//
//
// Original author KLG
//
// Notes:

// art includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"

// Mu2e includes.

#include "Mu2eG4/inc/constructCRV.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"
#include "Mu2eG4/inc/nestBox.hh"

// G4 includes

#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4Box.hh"

#include "Geant4/G4VSolid.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4PVPlacement.hh"

#include "Geant4/G4VisAttributes.hh"

#include "Geant4/G4RotationMatrix.hh"

#include "Geant4/G4SDManager.hh"

using namespace std;

namespace mu2e 
{
  void constructCRV( VolumeInfo const & parent, SimpleConfig const & _config)
  {
    GeomHandle<CosmicRayShield> CosmicRayShieldGeomHandle;

    Mu2eG4Helper& _helper       = *(art::ServiceHandle<Mu2eG4Helper>());
    AntiLeakRegistry& reg   = _helper.antiLeakRegistry();
    const auto& geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    
    geomOptions->loadEntry( _config, "crsVeto", "crs.veto" );
    const bool scintillatorShieldVisible   = geomOptions->isVisible("crsVeto"); 
    const bool scintillatorShieldDrawSolid = geomOptions->isSolid("crsVeto"); 
    const bool forceAuxEdgeVisible         = geomOptions->forceAuxEdgeVisible("crsVeto"); 
    const bool doSurfaceCheck              = geomOptions->doSurfaceCheck("crsVeto"); 
    //const bool placePV                     = geomOptions->placePV("crsVeto"); 
    
    int  verbosityLevel                    = _config.getInt("crs.verbosityLevel",0);
    bool forMARS                           = _config.getBool("crs.forMARS",false);

    std::string hallAirMaterialName = _config.getString("hall.insideMaterialName");

    CLHEP::Hep3Vector parentCenterInMu2e = parent.centerInMu2e();

    std::vector<CRSScintillatorShield> const &shields = CosmicRayShieldGeomHandle->getCRSScintillatorShields();
    std::vector<CRSScintillatorShield>::const_iterator ishield;
    for(ishield=shields.begin(); ishield!=shields.end(); ++ishield) 
    {
      CRSScintillatorShield const & shield = *ishield;
      std::string const & CRVsectorName = shield.getName();

      if(verbosityLevel > 0) cout << __func__ << " constructing            : " << CRVsectorName << endl;

      // all materials and dimensions are the same within a particular CRV sector
      CRSScintillatorBarDetail const & barDetail = shield.getCRSScintillatorBarDetail();

      G4Material* scintillatorBarMaterial = findMaterialOrThrow(barDetail.getMaterialName());
      std::vector<double> const &  scintillatorBarHalfLengths = barDetail.getHalfLengths();

      G4VSolid* scintillatorBarSolid = new G4Box(CRVsectorName,
                                                 scintillatorBarHalfLengths[0],
                                                 scintillatorBarHalfLengths[1],
                                                 scintillatorBarHalfLengths[2]);

      G4LogicalVolume* scintillatorBarLogical = new G4LogicalVolume(scintillatorBarSolid,
                                                                    scintillatorBarMaterial,
                                                                    CRVsectorName);

      // visibility attributes
      if (!scintillatorShieldVisible) 
      {
        scintillatorBarLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
      }
      else 
      {
        G4Colour  orange(.75, .55, .0);
        G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, orange));
        visAtt->SetForceSolid(scintillatorShieldDrawSolid);
        visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
        scintillatorBarLogical->SetVisAttributes(visAtt);
      }

      if (verbosityLevel > 1) 
      {
        G4SDManager::GetSDMpointer()->SetVerboseLevel(verbosityLevel-1);
      }

      // build one counter mother board, which will be the same for all counters of this CRV sector
      G4Material* CMBMaterial = findMaterialOrThrow(barDetail.getCMBMaterialName());
      std::vector<double> const & CMBHalfLengths = barDetail.getCMBHalfLengths();

      G4VSolid* CMBSolid = new G4Box("CRSCMB_"+CRVsectorName,
                                     CMBHalfLengths[0],
                                     CMBHalfLengths[1],
                                     CMBHalfLengths[2]);

      G4LogicalVolume* CMBLogical = new G4LogicalVolume(CMBSolid,CMBMaterial, "CRSCMB_"+CRVsectorName);

      // visibility attributes
      if (!scintillatorShieldVisible) 
      {
        CMBLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
      }
      else 
      {
        G4Colour  green(.0, .55, .55);
        G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, green));
        visAtt->SetForceSolid(scintillatorShieldDrawSolid);
        visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
        CMBLogical->SetVisAttributes(visAtt);
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
          const std::vector<double> &layerHalflengths=ilayer->getHalfLengths();

          //need the position relative to the parent (air) center
          const CLHEP::Hep3Vector &layerCenterInMu2e=ilayer->getPosition();
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
                                         false/*doSurfaceCheck*/);
	  if(doSurfaceCheck) checkForOverlaps(layerInfo.physical, _config, verbosityLevel>0);

          const std::vector<std::shared_ptr<CRSScintillatorBar> > &bars = ilayer->getBars();
          std::vector<std::shared_ptr<CRSScintillatorBar> >::const_iterator ibar;
          for(ibar=bars.begin(); ibar!=bars.end(); ++ibar) 
          {
            const CRSScintillatorBar &bar = **ibar; 

            //bar.getPosition() returns the bar position in Mu2e coordinates
            //need the position relative to the layer center
            CLHEP::Hep3Vector barLayerOffset = bar.getPosition() - layerCenterInMu2e; 
            CLHEP::Hep3Vector CMB0LayerOffset = bar.getCMBPosition(0) - layerCenterInMu2e; 
            CLHEP::Hep3Vector CMB1LayerOffset = bar.getCMBPosition(1) - layerCenterInMu2e; 

            G4VPhysicalVolume* pv = new G4PVPlacement(NULL,
                                                      barLayerOffset,
                                                      scintillatorBarLogical,
                                                      bar.name("CRSScintillatorBar_"),
                                                      layerInfo.logical,
                                                      false,
                                                      bar.index().asInt(),
                                                      false);
            if(doSurfaceCheck) checkForOverlaps(pv, _config, verbosityLevel>0);

            if(verbosityLevel > 4) std::cout << bar.name("CRSScintillatorBar_") << " " << bar.getPosition() <<std::endl;

//FIXME: this is temporary until the GDML issue is fixed
if(!_config.getBool("crs.hideCRVCMBs"))
{
            if(bar.hasCMB(0))
            {
              //MARS requires gdml file with unique logical volumes for the CMBs
              if(forMARS) CMBLogical = new G4LogicalVolume(CMBSolid,CMBMaterial, bar.name("CRSCMB0_")); 

              G4VPhysicalVolume* pvCMB0 = new G4PVPlacement(NULL,
                                                            CMB0LayerOffset,
                                                            CMBLogical,
                                                            bar.name("CRSCMB0_"),
                                                            layerInfo.logical,
                                                            false,
                                                            2*bar.index().asInt(),
                                                            false);
              if(doSurfaceCheck) checkForOverlaps(pvCMB0, _config, verbosityLevel>0);

              if(verbosityLevel > 4) std::cout << bar.name("CRSCMB0_") << " " << bar.getCMBPosition(0) <<std::endl;
            }

            if(bar.hasCMB(1))
            {
              //MARS requires gdml file with unique logical volumes for the CMBs
              if(forMARS) CMBLogical = new G4LogicalVolume(CMBSolid,CMBMaterial, bar.name("CRSCMB1_")); 

              G4VPhysicalVolume* pvCMB1 = new G4PVPlacement(NULL,
                                                            CMB1LayerOffset,
                                                            CMBLogical,
                                                            bar.name("CRSCMB1_"),
                                                            layerInfo.logical,
                                                            false,
                                                            2*bar.index().asInt()+1,
                                                            false);
              if(doSurfaceCheck) checkForOverlaps(pvCMB1, _config, verbosityLevel>0);

              if(verbosityLevel > 4) std::cout << bar.name("CRSCMB1_") << " " << bar.getCMBPosition(1) <<std::endl;
            }
}
          } //ibar
        } //ilayer

        //absorber sheets
        const std::vector<CRSAbsorberLayer> &absorberLayers = imodule->getAbsorberLayers();
        const std::string &absorberNameBase = imodule->name("CRSAbsorber_");
        for(unsigned int absorberLayerNumber=0; absorberLayerNumber<absorberLayers.size(); absorberLayerNumber++) 
        {
          const std::string absorberName = absorberNameBase+"_"+std::to_string(absorberLayerNumber);
          const std::vector<double> &absorberLayerHalflengths=absorberLayers[absorberLayerNumber].getHalfLengths();
          const CLHEP::Hep3Vector &absorberLayerCenterInMu2e=absorberLayers[absorberLayerNumber].getPosition();
          CLHEP::Hep3Vector absorberLayerAirOffset = absorberLayerCenterInMu2e - parentCenterInMu2e;

          const std::string &absorberMaterialName = shield.getAbsorberMaterialName();
          G4Material* absorberMaterial = findMaterialOrThrow(absorberMaterialName);

          G4VSolid* absorberSolid = new G4Box(absorberName,
                                              absorberLayerHalflengths[0],
                                              absorberLayerHalflengths[1],
                                              absorberLayerHalflengths[2]);

          G4LogicalVolume* absorberLogical = new G4LogicalVolume(absorberSolid,
                                                                 absorberMaterial,
                                                                 absorberName);

          if(!scintillatorShieldVisible) 
          {
            absorberLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
          }
          else 
          {
            G4Colour  darkorange  (.45, .25, .0);
            G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, darkorange));
            visAtt->SetForceSolid(scintillatorShieldDrawSolid);
            visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
            absorberLogical->SetVisAttributes(visAtt);
          }

          G4VPhysicalVolume* pv = new G4PVPlacement(NULL,
                                                    absorberLayerAirOffset,
                                                    absorberLogical,
                                                    absorberName,
                                                    parent.logical,
                                                    false,
                                                    0,
                                                    false);
          if(doSurfaceCheck) 
          {
            checkForOverlaps( pv, _config, verbosityLevel>0);
          }
        } //absorberLayerNumber


        //aluminum sheets
        const std::vector<CRSAluminumSheet> &aluminumSheets = imodule->getAluminumSheets();
        const std::string &aluminumSheetNameBase = imodule->name("CRSAluminumSheet_");
        for(unsigned int aluminumSheetNumber=0; aluminumSheetNumber<aluminumSheets.size(); aluminumSheetNumber++) 
        {
          const std::string aluminumSheetName = aluminumSheetNameBase+"_"+std::to_string(aluminumSheetNumber);
          const std::vector<double> &aluminumSheetHalflengths=aluminumSheets[aluminumSheetNumber].getHalfLengths();
          const CLHEP::Hep3Vector &aluminumSheetCenterInMu2e=aluminumSheets[aluminumSheetNumber].getPosition();
          CLHEP::Hep3Vector aluminumSheetAirOffset = aluminumSheetCenterInMu2e - parentCenterInMu2e;

          const std::string &aluminumSheetMaterialName = shield.getAluminumSheetMaterialName();
          G4Material* aluminumSheetMaterial = findMaterialOrThrow(aluminumSheetMaterialName);

          G4VSolid* aluminumSheetSolid = new G4Box(aluminumSheetName,
                                                   aluminumSheetHalflengths[0],
                                                   aluminumSheetHalflengths[1],
                                                   aluminumSheetHalflengths[2]);

          G4LogicalVolume* aluminumSheetLogical = new G4LogicalVolume(aluminumSheetSolid,
                                                                      aluminumSheetMaterial,
                                                                      aluminumSheetName);

          if(!scintillatorShieldVisible) 
          {
            aluminumSheetLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
          }
          else 
          {
            G4Colour  darkorange  (.45, .25, .0);
            G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, darkorange));
            visAtt->SetForceSolid(scintillatorShieldDrawSolid);
            visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
            aluminumSheetLogical->SetVisAttributes(visAtt);
          }

          G4VPhysicalVolume* pv = new G4PVPlacement(NULL,
                                                    aluminumSheetAirOffset,
                                                    aluminumSheetLogical,
                                                    aluminumSheetName,
                                                    parent.logical,
                                                    false,
                                                    0,
                                                    false);
          if(doSurfaceCheck) 
          {
            checkForOverlaps( pv, _config, verbosityLevel>0);
          }
        } //aluminumSheetNumber


        //FEBs
        const std::vector<CRSFEB> &FEBs = imodule->getFEBs();
        const std::string &FEBNameBase = imodule->name("CRSFEB_");
        for(unsigned int FEBNumber=0; FEBNumber<FEBs.size(); FEBNumber++) 
        {
          const std::string FEBName = FEBNameBase+"_"+std::to_string(FEBNumber);
          const std::vector<double> &FEBHalflengths=FEBs[FEBNumber].getHalfLengths();
          const CLHEP::Hep3Vector &FEBCenterInMu2e=FEBs[FEBNumber].getPosition();
          CLHEP::Hep3Vector FEBAirOffset = FEBCenterInMu2e - parentCenterInMu2e;

          const std::string &FEBMaterialName = shield.getFEBMaterialName();
          G4Material* FEBMaterial = findMaterialOrThrow(FEBMaterialName);

          G4VSolid* FEBSolid = new G4Box(FEBName,
                                         FEBHalflengths[0],
                                         FEBHalflengths[1],
                                         FEBHalflengths[2]);

          G4LogicalVolume* FEBLogical = new G4LogicalVolume(FEBSolid,
                                                            FEBMaterial,
                                                            FEBName);

          if(!scintillatorShieldVisible) 
          {
            FEBLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
          }
          else 
          {
            G4Colour  darkorange  (.45, .25, .0);
            G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, darkorange));
            visAtt->SetForceSolid(scintillatorShieldDrawSolid);
            visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
            FEBLogical->SetVisAttributes(visAtt);
          }

          G4VPhysicalVolume* pv = new G4PVPlacement(NULL,
                                                    FEBAirOffset,
                                                    FEBLogical,
                                                    FEBName,
                                                    parent.logical,
                                                    false,
                                                    0,
                                                    false);
          if(doSurfaceCheck) 
          {
            checkForOverlaps( pv, _config, verbosityLevel>0);
          }

          if(verbosityLevel > 4) std::cout << FEBName << " " << FEBCenterInMu2e <<std::endl;
        } //FEBNumber
      } //imodule
    } //ishield

    //steel support structures
    std::vector<CRSSupportStructure> const &supportStructures = CosmicRayShieldGeomHandle->getSupportStructures();
    std::vector<CRSSupportStructure>::const_iterator iSupportStructure;
    for(iSupportStructure=supportStructures.begin(); iSupportStructure!=supportStructures.end(); ++iSupportStructure) 
    {
      CRSSupportStructure const & supportStructure = *iSupportStructure;
      std::string const & name = supportStructure.getName();

      const std::vector<double> &halflengths=supportStructure.getHalfLengths();
      const CLHEP::Hep3Vector &positionInMu2e = supportStructure.getPosition();
      CLHEP::Hep3Vector supportStructureAirOffset = positionInMu2e - parentCenterInMu2e;

      const std::string &materialName = supportStructure.getMaterialName();
      G4Material* material = findMaterialOrThrow(materialName);

      G4VSolid* supportStructureSolid = new G4Box(name,
                                         halflengths[0],
                                         halflengths[1],
                                         halflengths[2]);

      G4LogicalVolume* supportStructureLogical = new G4LogicalVolume(supportStructureSolid, material, name);

      if(!scintillatorShieldVisible) 
      {
        supportStructureLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
      }
      else 
      {
        G4Colour  darkorange  (.45, .25, .0);
        G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, darkorange));
        visAtt->SetForceSolid(scintillatorShieldDrawSolid);
        visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
        supportStructureLogical->SetVisAttributes(visAtt);
      }

      G4VPhysicalVolume* pv = new G4PVPlacement(NULL,
                                                supportStructureAirOffset,
                                                supportStructureLogical,
                                                name,
                                                parent.logical,
                                                false,
                                                0,
                                                false);
      if(doSurfaceCheck) 
      {
        checkForOverlaps( pv, _config, verbosityLevel>0);
      }
    } //iSupportStructure

  } //construct CRV

} //namespace mu2e
