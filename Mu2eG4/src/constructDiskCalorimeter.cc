//
// Free function to create the Disk calorimeter 
//
//
// Original author Bertrand Echenard
//
//
// Note about steps in the disk. There are two ways to add the steps between the crystals and the disk edges. 
//      1) create each edge separately (by subtracting the inner hole or as the intersection between the outside edge and the bar)
//      2) make an union comprising of all the inner steps described as boxes, then subtract the hole inside. For the outer steps
//         describe the hole inside the disk as a sum of "empty" bars, then remove them from the cylinder.
//      
//   It turns out that version 1 is more annoying to code, but produced much faster code, so we picked it.  
//

#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "Mu2eG4/inc/constructDiskCalorimeter.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestTorus.hh"
#include "Mu2eG4/inc/CaloCrystalSD.hh"
#include "Mu2eG4/inc/CaloReadoutSD.hh"
#include "Mu2eG4/inc/CaloReadoutCardSD.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"

#include "G4Box.hh"
#include "G4Polyhedra.hh"
#include "G4Tubs.hh"
#include "G4Torus.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4VSolid.hh"
#include "G4String.hh"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/TwoVector.h"

#include <array>
#include <iostream>
#include <sstream>



namespace mu2e {

  VolumeInfo constructDiskCalorimeter(const VolumeInfo&  mother, const SimpleConfig& config,
                                      SensitiveDetectorHelper const& sdHelper)
  {

       const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
       geomOptions->loadEntry( config, "calorimeterOutline", "calorimeter.outline" );

       const bool isCalorimeterVisible = geomOptions->isVisible("calorimeterOutline"); 
       const bool isCalorimeterSolid   = geomOptions->isSolid("calorimeterOutline"); 
       const bool forceEdge            = config.getBool("g4.forceEdge",false);
       const bool doSurfaceCheck       = config.getBool("g4.doSurfaceCheck",false) || config.getBool("calorimeter.doSurfaceCheck",false);
       const int  verbosityLevel       = config.getInt("calorimeter.verbosityLevel",1);

       MaterialFinder materialFinder(config);
       G4Material* vacuumMaterial  = materialFinder.get("calorimeter.vacuumMaterial");

       G4PVPlacement* pv;


    //--------------------------------------
    // Construct calorimeter mother volume

       const DiskCalorimeter& cal  = *(GeomHandle<DiskCalorimeter>());
       
       const int  crateVersion   = config.getInt("calorimeter.crateVersion",1);
       const unsigned nDisks     = cal.nDisk();
       G4double mother_inRadius  = cal.caloInfo().getDouble("envelopeRadiusIn");
       G4double mother_outRadius = cal.caloInfo().getDouble("envelopeRadiusOut");
       G4double mother_z0        = cal.caloInfo().getDouble("envelopeZ0");
       G4double mother_z1        = cal.caloInfo().getDouble("envelopeZ1");
       G4double mother_zlength   = mother_z1-mother_z0;
       G4double mother_zCenter   = (mother_z1+mother_z0)/2.0;
       G4double vdMargin         = 0.5*CLHEP::mm;


       // Make the mother volume for the calorimeter.
       const CLHEP::Hep3Vector& posDS3  = mother.centerInMu2e();
       G4ThreeVector posCaloMother      = G4ThreeVector(posDS3.x(), 0, mother_zCenter);
       G4ThreeVector posCaloMotherInDS  = posCaloMother - posDS3;

       TubsParams caloParams(mother_inRadius,mother_outRadius,mother_zlength/2.0, 0., CLHEP::twopi);
       VolumeInfo calorimeterInfo = nestTubs("CalorimeterMother",caloParams,vacuumMaterial,0,posCaloMotherInDS,mother,0,
                                 	     isCalorimeterVisible,G4Colour::Blue(),isCalorimeterSolid,forceEdge,
                                	     true,doSurfaceCheck);

       if ( verbosityLevel > 0) 
       {
            double zhl = static_cast<G4Tubs*>(calorimeterInfo.solid)->GetZHalfLength();
	    double CalorimeterOffsetInMu2eZ = calorimeterInfo.centerInMu2e()[CLHEP::Hep3Vector::Z];
	    std::cout << __func__ << " Calorimeter mother center in Mu2e   : " << calorimeterInfo.centerInMu2e() << std::endl;
	    std::cout << __func__ << " Calorimeter mother Z extent in Mu2e : " << CalorimeterOffsetInMu2eZ - zhl << ", " << CalorimeterOffsetInMu2eZ + zhl << std::endl;
       }




    //--------------------------------------
    // Construct the full disk / FEB volumes, place the subcomponents inside and put them in the mother volumes
    // IMPORTANT: keep a 1mm buffer for virtual detectors, except for disk outer radius

       VolumeInfo calorimeterFEB[nDisks];
       VolumeInfo calorimeterDisk[nDisks];
              
       G4LogicalVolume* FEBLog = (crateVersion > 1) ? caloBuildFEB(config,materialFinder, cal ) : nullptr;
       G4Tubs* FEB             = (crateVersion > 1) ? static_cast<G4Tubs*>(FEBLog->GetSolid())  : nullptr;
               
       for (unsigned int idisk=0;idisk<nDisks;++idisk)
       {
           G4LogicalVolume* frontPlateLog = caloBuildFrontPlate(config,materialFinder, cal, idisk );
           G4LogicalVolume* diskLog       = caloBuildDisk(config,materialFinder, cal, idisk);
           G4LogicalVolume* backPlateLog  = caloBuildBackPlate(config,materialFinder, cal, idisk );

           G4Tubs* frontPlate = (frontPlateLog != nullptr) ? static_cast<G4Tubs*>(frontPlateLog->GetSolid()) : nullptr;
           G4Tubs* backPlate  = (backPlateLog  != nullptr) ? static_cast<G4Tubs*>(backPlateLog->GetSolid())  : nullptr; 
           G4Tubs* disk       = static_cast<G4Tubs*>(diskLog->GetSolid());

           G4double R0disk    = mother_inRadius +vdMargin;
           G4double R1disk    = mother_outRadius-vdMargin;
           G4double zHalfDisk = disk->GetZHalfLength();
           G4double zHalfFP   = (frontPlateLog != nullptr) ? frontPlate->GetZHalfLength() : 0;
           G4double zHalfBP   = (backPlateLog  != nullptr) ? backPlate->GetZHalfLength() : 0;
           G4double zHalftot  = zHalfFP+zHalfDisk+zHalfBP;
           
           double diskpar[5] = {R0disk,R1disk,zHalftot,0,CLHEP::twopi};
           std::ostringstream discname;  discname<<"caloDisk_" <<idisk;

	   //origin gives the position of the center of the disk, irrespective of the coordinate origin set in the calo description
	   G4ThreeVector posDisk = cal.disk(idisk).geomInfo().origin() - posCaloMother;

	   calorimeterDisk[idisk] = nestTubs(discname.str(),diskpar,vacuumMaterial,&cal.disk(idisk).geomInfo().rotation(),posDisk,
                                             calorimeterInfo,idisk,
                        		     isCalorimeterVisible,G4Colour::White(),0,forceEdge,true,1 );

           if (frontPlateLog != nullptr) pv = new G4PVPlacement(0,G4ThreeVector(0,0,-zHalftot+zHalfFP),frontPlateLog,"caloFP_PV",
                                                                calorimeterDisk[idisk].logical,false,0,false);
           if (frontPlateLog != nullptr) doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);

           pv = new G4PVPlacement(0,G4ThreeVector(0,0,-zHalftot+2*zHalfFP+zHalfDisk),diskLog,"caloCase_PV",
                                  calorimeterDisk[idisk].logical,false,0,false);
           doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);

           if (backPlateLog != nullptr) pv = new G4PVPlacement(0,G4ThreeVector(0,0,+zHalftot-zHalfBP),backPlateLog,"caloBP_PV",
                                                               calorimeterDisk[idisk].logical,false,0,false);
           if (backPlateLog != nullptr) doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);


	   if ( crateVersion > 1 )
           {
   	      std::ostringstream cratename; cratename<<"caloFEB_" <<idisk;
              double FEBpar[5] = {FEB->GetInnerRadius()-vdMargin,FEB->GetOuterRadius()+vdMargin,FEB->GetZHalfLength()+0.5*vdMargin,
                                  FEB->GetStartPhiAngle(),FEB->GetDeltaPhiAngle()};

              G4ThreeVector posFEB  = posCaloMotherInDS + cal.disk(idisk).geomInfo().originLocal() + G4ThreeVector(0,0,cal.disk(idisk).geomInfo().crateDeltaZ());

	      calorimeterFEB[idisk] = nestTubs(cratename.str(),FEBpar,vacuumMaterial,&cal.disk(idisk).geomInfo().rotation(),posFEB,
                        		       mother,idisk,
                                               isCalorimeterVisible,G4Colour::White(),0,forceEdge,true,doSurfaceCheck );

              pv = new G4PVPlacement(0,G4ThreeVector(0,0,0),FEBLog,"caloFEB_PV",calorimeterFEB[idisk].logical,false,0,false);
              doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
           }

	   if ( verbosityLevel > 0) 
	   {
	       std::cout << __func__ << " CalorimeterDisk center in Mu2e    : " << calorimeterDisk[idisk].centerInMu2e() << std::endl;
	       std::cout << __func__ << " CalorimeterDisk Z extent in Mu2e  : " << calorimeterDisk[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] - zHalftot << ", " << calorimeterDisk[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] + zHalftot << std::endl;
	       std::cout << __func__ << " Calorimeter FP / DISK / BP half depth: "<<zHalfFP<<" / "<<zHalfDisk<<" / "<<zHalfBP<<std::endl;
           }

	}

        if ( verbosityLevel > 0) std::cout << __func__ << " Calorimeter constructed "<<std::endl;
        return calorimeterInfo;
  }
  




  //--------------------------------------------------------------------------------------------------------------------------------
  // Front plate: a sandwich of carbon fiberplate - foam plane - carbon fiber. The pipes are inside the foam layer, closest to the beam. 
  // the big pipe is aligned with the small pipes (centered on the same z). Big cooling pipe needs to be larger than 
  // carbon thickness+small pipe radius and smaller than Carbon thick + foam thick - small pipe radius or the model needs to be updated.
  //  
  G4LogicalVolume* caloBuildFrontPlate(const SimpleConfig& config, MaterialFinder& materialFinder, const DiskCalorimeter& cal, int idisk)
  { 
       const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
       geomOptions->loadEntry( config, "calorimeterPipe", "calorimeter.pipe" );

       const bool isPipeVisible          = geomOptions->isVisible("calorimeterPipe"); 
       const bool isPipeSolid            = geomOptions->isSolid("calorimeterPipe"); 
       const bool forceEdge              = config.getBool("g4.forceEdge",false);
       const bool doSurfaceCheck         = config.getBool("g4.doSurfaceCheck",false) || config.getBool("calorimeter.doSurfaceCheck",false);
       const int  verbosityLevel         = config.getInt("calorimeter.verbosityLevel",1);
       G4VPhysicalVolume* pv;

       G4Material* vacuumMaterial        = materialFinder.get("calorimeter.vacuumMaterial");
       G4Material* FPFoamMaterial        = materialFinder.get("calorimeter.FPFoamMaterial");
       G4Material* FPCarbonMaterial      = materialFinder.get("calorimeter.FPCarbonMaterial");
       G4Material* coolPipeMaterial      = materialFinder.get("calorimeter.coolPipeMaterial");
       G4Material* pipeMaterial          = materialFinder.get("calorimeter.pipeMaterial");

       G4double FPInnerRadius            = cal.caloInfo().getDouble("FPInnerRadius");
       G4double FPOuterRadius            = cal.caloInfo().getDouble("FPOuterRadius");
       G4double FPCarbonHalfThick        = cal.caloInfo().getDouble("FPCarbonZLength")/2.0;  
       G4double FPFoamHalfThick          = cal.caloInfo().getDouble("FPFoamZLength")/2.0;  
       G4double FPCoolPipeTorRadius      = cal.caloInfo().getDouble("FPCoolPipeTorRadius");  
       G4double FPCoolPipeRadius         = cal.caloInfo().getDouble("FPCoolPipeRadius");  
       G4double FPCoolPipeThickness      = cal.caloInfo().getDouble("FPCoolPipeThickness");

       G4int nPipes                      = cal.caloInfo().getInt("nPipes");      
       G4double pipeRadius               = cal.caloInfo().getDouble("pipeRadius");
       G4double pipeThickness            = cal.caloInfo().getDouble("pipeThickness");
       G4double pipeInitSeparation       = cal.caloInfo().getDouble("pipeInitSeparation");    
       std::vector<double> pipeTorRadius = cal.caloInfo().getVDouble("pipeTorRadius");
              
       G4double frontPanelHalfThick      = (2.0*FPCarbonHalfThick+2.0*FPFoamHalfThick-pipeRadius+FPCoolPipeRadius)/2.0;
       G4double ZposCarbon2              = frontPanelHalfThick-FPCarbonHalfThick;
       G4double ZposFoam                 = ZposCarbon2-FPCarbonHalfThick-FPFoamHalfThick;
       G4double ZposCarbon1              = ZposFoam-FPFoamHalfThick-FPCarbonHalfThick;
       G4double ZposPipe                 = frontPanelHalfThick-2*FPCarbonHalfThick-2*FPFoamHalfThick+pipeRadius;
      
       if (nPipes==0) return nullptr;

       //this is the full front panel
       G4Tubs* frontPlate = new G4Tubs("caloFrontPlate",FPInnerRadius,FPCoolPipeTorRadius+FPCoolPipeRadius,frontPanelHalfThick,0,CLHEP::twopi);       
       G4LogicalVolume* frontPlateLog = caloBuildLogical(frontPlate, vacuumMaterial, "caloFrontPlateLog",0,G4Color::White(),0,0);

       //carbon fiber panels
       G4Tubs*          frontPanelCarb    = new G4Tubs("caloFPCarb",FPInnerRadius,FPOuterRadius,FPCarbonHalfThick,0,CLHEP::twopi);       
       G4LogicalVolume* frontPanelCarbLog = caloBuildLogical(frontPanelCarb, FPCarbonMaterial, "caloFPCarbLog",isPipeVisible,G4Color::Grey(),0,forceEdge);       
       
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,ZposCarbon1), frontPanelCarbLog, "caloFPCarbPV1", frontPlateLog, false, 0, false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);                
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,ZposCarbon2), frontPanelCarbLog, "caloFPCarbPV2", frontPlateLog, false, 0, false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);                

       //Foam panel
       G4Tubs*          frontPanelFoam    = new G4Tubs("caloFPFoam",FPInnerRadius,FPOuterRadius,FPFoamHalfThick,0,CLHEP::twopi);       
       G4LogicalVolume* frontPanelFoamLog = caloBuildLogical(frontPanelFoam, FPFoamMaterial, "caloFPFoamLog",isPipeVisible,G4Color::Brown(),0,forceEdge);
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,ZposFoam), frontPanelFoamLog, "caloFPFoamPV", frontPlateLog, false, 0, false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);                
                    
       //cooling pipes on the edge
       double pipeTotalSep = pipeTorRadius.back()-pipeTorRadius.front();
       double angMax = CLHEP::pi/2.0-std::asin((pipeInitSeparation+pipeTotalSep)/FPOuterRadius)-0.1;
       G4RotationMatrix* rotFPPipe = new G4RotationMatrix(CLHEP::HepRotation::IDENTITY);
       rotFPPipe->rotateY(CLHEP::pi);
       
       G4Torus*         coolFP     = new G4Torus("caloCoolFP",FPCoolPipeRadius-FPCoolPipeThickness, FPCoolPipeRadius, FPCoolPipeTorRadius, angMax, CLHEP::twopi-2.0*angMax);
       G4LogicalVolume* coolFPLog  = caloBuildLogical(coolFP, coolPipeMaterial, "caloCoolFPLog",isPipeVisible,G4Color::Red(),isPipeSolid,0);
       pv = new G4PVPlacement(rotFPPipe,G4ThreeVector(0.0,0.0,ZposPipe), coolFPLog, "caloCoolFPPV", frontPlateLog, false, 0, false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);                
       

       //pipes nside foam
       G4RotationMatrix* rotPipe1 = new G4RotationMatrix(CLHEP::HepRotation::IDENTITY);
       rotPipe1->rotateZ(CLHEP::pi/2.0);
       G4RotationMatrix* rotPipe2 = new G4RotationMatrix(CLHEP::HepRotation::IDENTITY);
       rotPipe2->rotateZ(1.5*CLHEP::pi);
       G4RotationMatrix* rotPipeFlat = new G4RotationMatrix(CLHEP::HepRotation::IDENTITY);
       rotPipeFlat->rotateX(CLHEP::pi/2.0);

       for (int ipipe=0; ipipe<nPipes; ++ipipe)
       {
           double xpipe  = pipeInitSeparation+pipeTorRadius[ipipe]-pipeTorRadius[0];
           double angle  = std::asin(xpipe/pipeTorRadius[ipipe]); //angle taken w.r.t y axis!
           double length = sqrt(FPOuterRadius*FPOuterRadius-xpipe*xpipe) - (pipeTorRadius[ipipe]+pipeRadius)*cos(angle) - 2.0*pipeRadius;
           double y0     = (pipeTorRadius[ipipe]+pipeRadius)*cos(angle)+0.5*length;
           double y1     = -(pipeTorRadius[ipipe]+pipeRadius)*cos(angle)-0.5*length;
           double z      = -FPFoamHalfThick+pipeRadius;

           G4Torus* pipe1 = new G4Torus("caloPipe",pipeRadius-pipeThickness, pipeRadius, pipeTorRadius[ipipe],angle, CLHEP::pi-2*angle);
           G4Tubs*  pipe2 = new G4Tubs("caloPipe2",pipeRadius-pipeThickness, pipeRadius,0.5*length,0,CLHEP::twopi);       
           G4LogicalVolume* pipe1Log = caloBuildLogical(pipe1, pipeMaterial, "caloPipe1Log",isPipeVisible,G4Color::Cyan(),isPipeSolid,forceEdge);
           G4LogicalVolume* pipe2Log = caloBuildLogical(pipe2, pipeMaterial, "caloPipe2Log",isPipeVisible,G4Color::Cyan(),isPipeSolid,forceEdge);

           pv = new G4PVPlacement(rotPipe1,G4ThreeVector(0,0,z), pipe1Log, "caloPipePV", frontPanelFoamLog, false, ipipe, false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);                
           pv = new G4PVPlacement(rotPipe2,G4ThreeVector(0,0,z), pipe1Log, "caloPipePV", frontPanelFoamLog, false, ipipe, false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);                                
           pv = new G4PVPlacement(rotPipeFlat,G4ThreeVector(xpipe, y0,z), pipe2Log, "caloPipePV", frontPanelFoamLog, false, ipipe , false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);           
           pv = new G4PVPlacement(rotPipeFlat,G4ThreeVector(xpipe,-y0,z), pipe2Log, "caloPipePV", frontPanelFoamLog, false, ipipe , false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);           
           pv = new G4PVPlacement(rotPipeFlat,G4ThreeVector(-xpipe,y1,z), pipe2Log, "caloPipePV", frontPanelFoamLog, false, ipipe , false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);           
           pv = new G4PVPlacement(rotPipeFlat,G4ThreeVector(-xpipe,-y1,z), pipe2Log, "caloPipePV", frontPanelFoamLog, false, ipipe , false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);           
       }

       return frontPlateLog;
   }




  //--------------------------------------------------------------------------------------------------------------------------------
  //construct central part of the disk with crystals
  G4LogicalVolume* caloBuildDisk(const SimpleConfig& config, MaterialFinder& materialFinder, const DiskCalorimeter& cal, int idisk)
  {
       const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
       geomOptions->loadEntry( config, "calorimeterCase", "calorimeter.case" );
       geomOptions->loadEntry( config, "calorimeterCrystal", "calorimeter.crystal" );

       const bool isDiskVisible          = geomOptions->isVisible("calorimeterCase"); 
       const bool isDiskSolid            = geomOptions->isSolid("calorimeterCase"); 
       const bool isCrystalVisible       = geomOptions->isVisible("calorimeterCrystal"); 
       const bool isCrystalSolid         = geomOptions->isSolid("calorimeterCrystal");        
       const bool forceEdge              = config.getBool("g4.forceEdge",false);
       const bool doSurfaceCheck         = config.getBool("g4.doSurfaceCheck",false) || config.getBool("calorimeter.doSurfaceCheck",false);
       const int  verbosityLevel         = config.getInt("calorimeter.verbosityLevel",1);
       G4VPhysicalVolume* pv;

       G4Material* vacuumMaterial        = materialFinder.get("calorimeter.vacuumMaterial");
       G4Material* crysMaterial          = materialFinder.get("calorimeter.crystalMaterial");
       G4Material* crysFrameMaterial     = materialFinder.get("calorimeter.crystalFrameMaterial");
       G4Material* wrapMaterial          = materialFinder.get("calorimeter.wrapperMaterial");    
       G4Material* innerRingMaterial     = materialFinder.get("calorimeter.innerRingMaterial");
       G4Material* outerRingMaterial     = materialFinder.get("calorimeter.outerRingMaterial");
       G4Material* coolPipeMaterial      = materialFinder.get("calorimeter.coolPipeMaterial");


       G4double crystalHalfXY            = cal.caloInfo().getDouble("crystalXYLength")/2.0;
       G4double crystalHalfZ             = cal.caloInfo().getDouble("crystalZLength")/2.0;    
       G4double crystalFrameHalfZ        = cal.caloInfo().getDouble("crystalFrameZLength")/2.0;    
       G4double crystalFrameThick        = cal.caloInfo().getDouble("crystalFrameThickness");    
       G4double wrapperHalfThick         = cal.caloInfo().getDouble("wrapperThickness")/2.0;    
       G4double wrapperHalfXY            = crystalHalfXY + 2*wrapperHalfThick;
       G4double wrapperHalfZ             = crystalHalfZ+2.0*crystalFrameHalfZ;
       G4double crystalFrameHalfXY       = crystalHalfXY-crystalFrameThick;

       G4double diskCaseHalfZLength      = cal.caloInfo().getDouble("diskCaseZLength")/2.0;
       G4double diskCaseInnerRadius      = cal.caloInfo().getDouble("diskCaseRadiusIn");
       G4double diskCaseOuterRadius      = cal.caloInfo().getDouble("diskCaseRadiusOut");
       G4double innerRingThickness       = cal.caloInfo().getDouble("diskInnerRingThickness");
       G4double outerRingThickness       = cal.caloInfo().getDouble("diskOuterRingThickness");
       G4double outerRingEdgeHalfZLength = cal.caloInfo().getDouble("diskOutRingEdgeZLength")/2.0;
       G4double outerRingEdgeThickness   = cal.caloInfo().getDouble("diskOutRingEdgeRLength");     
       G4double FPCoolPipeRadius         = cal.caloInfo().getDouble("FPCoolPipeRadius");  
       G4double FPCoolPipeThickness      = cal.caloInfo().getDouble("FPCoolPipeThickness");
       G4double coolPipeZpos             = (diskCaseHalfZLength - 2*FPCoolPipeRadius - 2.0*outerRingEdgeHalfZLength)/2.0 + FPCoolPipeRadius;

       std::vector<double> stepsInX      = cal.caloInfo().getVDouble("stepsInsideX");
       std::vector<double> stepsInY      = cal.caloInfo().getVDouble("stepsInsideY");
       std::vector<double> stepsOutX     = cal.caloInfo().getVDouble("stepsOutsideX");
       std::vector<double> stepsOutY     = cal.caloInfo().getVDouble("stepsOutsideY");
       G4double diskStepThickness        = cal.caloInfo().getDouble("diskStepThickness");     

       //-----------------------------------------
       // Build wrapper crystal unit, including frame. the wrapper is not covering the front face
       //----
       G4Box* crystalWrap    = new G4Box("caloCrystalWrap",   wrapperHalfXY,wrapperHalfXY,wrapperHalfZ);
       G4Box* crystal        = new G4Box("caloCrystal",       crystalHalfXY,crystalHalfXY,crystalHalfZ);
       G4Box* crystalFrame   = new G4Box("caloCrystalFrame",  crystalHalfXY,crystalHalfXY,crystalFrameHalfZ);
       G4Box* crystalFrameIn = new G4Box("caloCrystalFrameIn",crystalFrameHalfXY,crystalFrameHalfXY,crystalFrameHalfZ);
       G4SubtractionSolid* crystalHollowFrame = new G4SubtractionSolid("caloCrystalFullFrame", crystalFrame, crystalFrameIn, 0, G4ThreeVector(0,0,0));           

       G4LogicalVolume *crystalLog        = caloBuildLogical(crystal,            crysMaterial,      "caloCrystalLog",0, G4Color::Cyan(),0,0);
       G4LogicalVolume *crystalFrameLog   = caloBuildLogical(crystalHollowFrame, crysFrameMaterial, "caloCrystalFrameLog",isCrystalVisible, G4Color::White(),isCrystalSolid,0);
       G4LogicalVolume *wrapperLog        = caloBuildLogical(crystalWrap,        wrapMaterial,      "caloCrystalWrapLog",isCrystalVisible,G4Color::Cyan(),isCrystalSolid,forceEdge);  

       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,0.0),crystalLog,"caloCrysPV",wrapperLog,false,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,-wrapperHalfZ+crystalFrameHalfZ),crystalFrameLog,"caloCrysFrame1PV",wrapperLog,false,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,wrapperHalfZ-crystalFrameHalfZ),crystalFrameLog,"caloCrysFrame2PV",wrapperLog,false,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

       G4VSensitiveDetector* ccSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloCrystal());
       if (ccSD) crystalLog->SetSensitiveDetector(ccSD);


       //------------------------------------------------------------
       // Build disk inner ring, case ring (containing crystals) and outer ring
       G4Tubs* fullCrystalDisk  = new G4Tubs("calofullCrystalDisk", diskCaseInnerRadius,diskCaseOuterRadius+outerRingEdgeThickness,diskCaseHalfZLength,      0,CLHEP::twopi);       
       G4Tubs* innerRingDisk    = new G4Tubs("caloInnerRingDisk",   diskCaseInnerRadius,diskCaseInnerRadius+innerRingThickness,    diskCaseHalfZLength,      0,CLHEP::twopi);              
       G4Tubs* outerRingDisk    = new G4Tubs("caloOuterRingDisk",   diskCaseOuterRadius-outerRingThickness,diskCaseOuterRadius,    diskCaseHalfZLength,      0,CLHEP::twopi);       
       G4Tubs* outerRailDisk    = new G4Tubs("caloOuterRailCase",   diskCaseOuterRadius,diskCaseOuterRadius+outerRingEdgeThickness,outerRingEdgeHalfZLength, 0,CLHEP::twopi);       

       G4LogicalVolume* fullCrystalDiskLog = caloBuildLogical(fullCrystalDisk, vacuumMaterial,    "calofullCrystalDiskLog", 0, G4Color::Black(),0,0);   
       G4LogicalVolume* innerRingDiskLog   = caloBuildLogical(innerRingDisk,   innerRingMaterial, "caloInnerRingLog"      , isDiskVisible,G4Color::Yellow(),isDiskSolid,forceEdge);   
       G4LogicalVolume* outerRingDiskLog   = caloBuildLogical(outerRingDisk,   outerRingMaterial, "caloOuterRingLog"      , isDiskVisible,G4Color::Yellow(),isDiskSolid,forceEdge);   
       G4LogicalVolume* outerRailDiskLog   = caloBuildLogical(outerRailDisk,   outerRingMaterial, "caloouterRailDiskLog"  , isDiskVisible,G4Color::Yellow(),isDiskSolid,forceEdge);   

       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0,0),innerRingDiskLog,"innerRingDiskPV",fullCrystalDiskLog,true,0,false);
       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0,0),outerRingDiskLog,"outerRingDiskPV",fullCrystalDiskLog,true,0,false);
       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0,diskCaseHalfZLength-outerRingEdgeHalfZLength),outerRailDiskLog,"caloDiskRingCase1PV",fullCrystalDiskLog,true,0,false);
       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0,-diskCaseHalfZLength+outerRingEdgeHalfZLength),outerRailDiskLog,"caloDiskRingCase2PV",fullCrystalDiskLog,true,0,false);
       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
       
       
       //------------------------------------------------------------
       // add two cooling pipes running between the rails. Make them similar to front plate cooling pipes.        
       G4Torus* coolPipe            = new G4Torus("caloPipe",FPCoolPipeRadius-FPCoolPipeThickness, FPCoolPipeRadius, diskCaseOuterRadius+FPCoolPipeRadius,0,1.2*CLHEP::pi);
       G4LogicalVolume* coolPipeLog = caloBuildLogical(coolPipe, coolPipeMaterial, "caloCoolFPLog",isDiskVisible,G4Color::Red(),isDiskSolid,0);
       G4RotationMatrix* rotPipe    = new G4RotationMatrix(CLHEP::HepRotation::IDENTITY); rotPipe->rotateZ(0.1*CLHEP::pi);

       pv = new G4PVPlacement(rotPipe,G4ThreeVector(0.0,0.0,coolPipeZpos),  coolPipeLog, "caloCoolDisk1PV", fullCrystalDiskLog, false, 0, false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);                
       pv = new G4PVPlacement(rotPipe,G4ThreeVector(0.0,0.0,-coolPipeZpos), coolPipeLog, "caloCoolDisk2PV", fullCrystalDiskLog, false, 0, false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);                
       
       //--------------------------------------------------------------------
       // add the steps in the inner / outer part of the disk (see top Note). This is surprisingly compact...
       G4Tubs* holeIn = new G4Tubs("InnerCrystalEdge",0,diskCaseInnerRadius+innerRingThickness,diskCaseHalfZLength,0,CLHEP::twopi);           
       for (unsigned i=0;i<stepsInX.size();++i)
       {            
           G4Box* box  = new G4Box("ShimIn",stepsInX[i],wrapperHalfXY,crystalHalfZ);
           G4Box* box2 = new G4Box("ShimIn",stepsInX[i]-diskStepThickness,wrapperHalfXY-diskStepThickness,crystalHalfZ-diskStepThickness);
           G4SubtractionSolid* hollowBox = new G4SubtractionSolid("HollowBox shim", box, box2, 0, G4ThreeVector(0,0,0));           
           G4SubtractionSolid* cutShim = new G4SubtractionSolid("Cropped shim", hollowBox, holeIn,0, G4ThreeVector(0,-stepsInY[i],0));
           G4LogicalVolume* cutShimLog = caloBuildLogical(cutShim, innerRingMaterial, "cutShimLog",  isDiskVisible, G4Color::Yellow(), isDiskSolid, forceEdge);   
           pv = new G4PVPlacement(0,CLHEP::Hep3Vector(0,stepsInY[i],0),cutShimLog,"steps",fullCrystalDiskLog,true,0,false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       }

       G4Tubs* crystalOut = new G4Tubs("OuterCrystalEdge",0,diskCaseOuterRadius-outerRingThickness,diskCaseHalfZLength,0,CLHEP::twopi);           
       
       for (unsigned i=0;i<stepsOutX.size();++i)
       {
           double xlen = sqrt(diskCaseOuterRadius*diskCaseOuterRadius-stepsOutY[i]*stepsOutY[i])-stepsOutX[i];
           G4Box* box = new G4Box("ShimOut",xlen,wrapperHalfXY,crystalHalfZ);
           G4IntersectionSolid* cutShim1 = new G4IntersectionSolid("Cropped shim", crystalOut, box,0, G4ThreeVector(stepsOutX[i]+xlen,stepsOutY[i],0));
           G4LogicalVolume* cutShim1Log = caloBuildLogical(cutShim1, outerRingMaterial, "cutShimLog",  isDiskVisible, G4Color::Yellow(), isDiskSolid, forceEdge);   
           pv = new G4PVPlacement(0,CLHEP::Hep3Vector(0,0,0),cutShim1Log,"steps",fullCrystalDiskLog,true,0,false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

           G4IntersectionSolid* cutShim2 = new G4IntersectionSolid("Cropped shim", crystalOut, box,0, G4ThreeVector(-stepsOutX[i]-xlen,stepsOutY[i],0));
           G4LogicalVolume* cutShim2Log = caloBuildLogical(cutShim2, outerRingMaterial, "cutShimLog",  isDiskVisible, G4Color::Yellow(), isDiskSolid, forceEdge);   
           pv = new G4PVPlacement(0,CLHEP::Hep3Vector(0,0,0),cutShim2Log,"steps",fullCrystalDiskLog,true,0,false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       }

       //add the two pieces at the top and bottom       
       if (stepsOutY.size()>0)
       {
           double ylow = *std::min_element(stepsOutY.begin(),stepsOutY.end())-3.0*wrapperHalfXY;
           double yup  = *std::max_element(stepsOutY.begin(),stepsOutY.end())+3.0*wrapperHalfXY;

           G4Box* box = new G4Box("ShimOut",0.4*diskCaseOuterRadius,2.0*wrapperHalfXY,crystalHalfZ);
           G4IntersectionSolid* cutShim1 = new G4IntersectionSolid("Cropped shim", crystalOut, box,0, G4ThreeVector(0,ylow,0));
           G4IntersectionSolid* cutShim2 = new G4IntersectionSolid("Cropped shim", crystalOut, box,0, G4ThreeVector(0,yup,0));

           G4LogicalVolume* cutShim1Log = caloBuildLogical(cutShim1, outerRingMaterial, "cutShimLog",  isDiskVisible, G4Color::Yellow(), isDiskSolid, forceEdge);   
           G4LogicalVolume* cutShim2Log = caloBuildLogical(cutShim2, outerRingMaterial, "cutShimLog",  isDiskVisible, G4Color::Yellow(), isDiskSolid, forceEdge);   

           pv = new G4PVPlacement(0,CLHEP::Hep3Vector(0,0,0),cutShim1Log,"steps",fullCrystalDiskLog,true,0,false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
           pv = new G4PVPlacement(0,CLHEP::Hep3Vector(0,0,0),cutShim2Log,"steps",fullCrystalDiskLog,true,0,false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       }
      
       
       
       //------------------------------------------------------------
       // finally put the crystals inside the crystalDiskCase 
       int nTotCrystal(0);
       for (int i=0;i<idisk;++i) nTotCrystal+=cal.disk(idisk).nCrystals();
       
       for (unsigned ic=0; ic <cal.disk(idisk).nCrystals(); ++ic)
       {	      
	   G4int id = nTotCrystal+ic;
	   std::ostringstream name;name<<"caloCrystalPV_" <<id;
           CLHEP::Hep3Vector crystalPosition = cal.disk(idisk).crystal(ic).localPosition();
           crystalPosition.setZ(diskCaseHalfZLength-wrapperHalfZ);

           //if needed, move the crystal construction here for individual dimensions             

           pv = new G4PVPlacement(0,crystalPosition,wrapperLog,name.str(),fullCrystalDiskLog,true,id,false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       }
              
       return fullCrystalDiskLog;
   }
  
  
  
  
  //--------------------------------------------------------------------------------------------------------------------------------
  // build full backplate - yes this was annoying
  G4LogicalVolume* caloBuildBackPlate(const SimpleConfig& config, MaterialFinder& materialFinder, const DiskCalorimeter& cal, int idisk)
  {            
       const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
       geomOptions->loadEntry( config, "calorimeterRO", "calorimeter.readout" );

       const bool isROVisible           = geomOptions->isVisible("calorimeterRO"); 
       const bool isROSolid             = geomOptions->isSolid("calorimeterRO"); 
       const bool forceEdge             = config.getBool("g4.forceEdge",false);
       const bool doSurfaceCheck        = config.getBool("g4.doSurfaceCheck",false) || config.getBool("calorimeter.doSurfaceCheck",false);
       const int  verbosityLevel        = config.getInt("calorimeter.verbosityLevel",1);

       G4Material* vacuumMaterial       = materialFinder.get("calorimeter.vacuumMaterial");
       G4Material* ROMaterial           = materialFinder.get("calorimeter.readoutMaterial");
       G4Material* FEEMaterial          = materialFinder.get("calorimeter.FEEMaterial");
       G4Material* FEEBoxMaterial       = materialFinder.get("calorimeter.FEEBoxMaterial");
       G4Material* backPlateMaterial    = materialFinder.get("calorimeter.BackPlateMaterial");
       G4Material* pipeMaterial         = materialFinder.get("calorimeter.coolPipeMaterial");
       G4Material* stripMaterial        = materialFinder.get("calorimeter.BPStripMaterial");

       G4double BPInnerRadius           = cal.caloInfo().getDouble("FPInnerRadius"); //same as front plate
       G4double BPOuterRadius           = cal.caloInfo().getDouble("BPOuterRadius"); 
      
       G4double crystalHalfXY           = cal.caloInfo().getDouble("crystalXYLength")/2.0;
       G4double wrapperHalfXY           = crystalHalfXY + cal.caloInfo().getDouble("wrapperThickness");
       G4double ROHalfX                 = cal.caloInfo().getDouble("readoutXLength")/2.0;
       G4double ROHalfY                 = cal.caloInfo().getDouble("readoutYLength")/2.0;
       G4double ROHalfZ                 = cal.caloInfo().getDouble("readoutZLength")/2.0;  
       G4double holeHalfX               = cal.caloInfo().getDouble("BPHoleXLength")/2.0;
       G4double holeHalfY               = cal.caloInfo().getDouble("BPHoleYLength")/2.0;
       G4double holeHalfZ               = cal.caloInfo().getDouble("BPHoleZLength")/2.0;
       G4double stripHalfY              = wrapperHalfXY-holeHalfY-1;
       G4double stripHalfZ              = cal.caloInfo().getDouble("BPStripThickness")/2.0;
       G4double FEEHalfX                = cal.caloInfo().getDouble("FEEXLength")/2.0;
       G4double FEEHalfY                = cal.caloInfo().getDouble("FEEYLength")/2.0;
       G4double FEEHalfZ                = cal.caloInfo().getDouble("FEEZLength")/2.0;

       G4double FEEBoxThickness         = cal.caloInfo().getDouble("FEEBoxThickness");
       G4double FEEBoxHalfX             = holeHalfX + 2*FEEBoxThickness;
       G4double FEEBoxHalfY             = FEEHalfY + 2*FEEBoxThickness;
       G4double FEEBoxHalfZ             = FEEHalfZ + 2*FEEBoxThickness;
       
       G4double BPPipeRadiusHigh        = cal.caloInfo().getDouble("BPPipeRadiusHigh");
       G4double BPPipeRadiusLow         = cal.caloInfo().getDouble("BPPipeRadiusLow");
       G4double BPPipeThickness         = cal.caloInfo().getDouble("BPPipeThickness");
       G4double BPPipeTorRadiusHigh     = BPOuterRadius - BPPipeRadiusHigh; 
       G4double BPPipeTorRadiusLow      = BPOuterRadius - 3.0*BPPipeRadiusHigh; 
       G4double BPPipeHalfZOffset       = cal.caloInfo().getDouble("BPPipeZOffset")/2.0;
       G4double BPFEEHalfZ              = FEEBoxHalfZ + BPPipeHalfZOffset + BPPipeRadiusHigh;

       std::vector<double> stepsInX   = cal.caloInfo().getVDouble("stepsInsideX");
       std::vector<double> stepsInY   = cal.caloInfo().getVDouble("stepsInsideY");
       std::vector<double> stepsOutX  = cal.caloInfo().getVDouble("stepsOutsideX");
       std::vector<double> stepsOutY  = cal.caloInfo().getVDouble("stepsOutsideY");
       int nstepsInX                  = int(stepsInX.size());
       int nstepsOutX                 = int(stepsOutX.size());
 

       G4VPhysicalVolume* pv;

       if (cal.caloInfo().nROPerCrystal()==0) return nullptr;


       //-------------------------------------------------------------------------------------------
       // Build hole in back plate, place SiPM at bottom of hole
       //----
       G4Box *holeBack    = new G4Box("caloHoleBack",    holeHalfX,holeHalfY,holeHalfZ);
       G4Box *crystalRO   = new G4Box("caloCrystalRO",   ROHalfX,ROHalfY,ROHalfZ);

       G4LogicalVolume* holeBackLog    = caloBuildLogical(holeBack,   vacuumMaterial, "caloHoleBackLog",  isROVisible, G4Color::Black(),isROSolid,forceEdge);   
       G4LogicalVolume* crystalROLog   = caloBuildLogical(crystalRO,  ROMaterial,     "caloCrystalROLog", isROVisible, G4Color::Grey(), isROSolid,forceEdge);        

       pv = new G4PVPlacement(0,G4ThreeVector(ROHalfX, 0, -holeHalfZ+ROHalfZ),  crystalROLog,"caloROPV_0", holeBackLog, true,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(-ROHalfX, 0, -holeHalfZ+ROHalfZ), crystalROLog,"caloROPV_1", holeBackLog, true,1,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

       // add sensitive detector    
       G4VSensitiveDetector* crSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloReadout());
       if (crSD) crystalROLog->SetSensitiveDetector(crSD);


       //----------------------
       // Build FEE in copper box
       //----
       G4Box *FEEBox    = new G4Box("caloFEEBox",  FEEBoxHalfX,FEEBoxHalfY,FEEBoxHalfZ);
       G4Box *FEEBoxIn  = new G4Box("caloFEEBoxIn",FEEBoxHalfX-FEEBoxThickness,FEEBoxHalfY-FEEBoxThickness,FEEBoxHalfZ-0.5*FEEBoxThickness);
       G4Box *FEECard   = new G4Box("caloFEECard", FEEHalfX,FEEHalfY,FEEHalfZ);

       G4LogicalVolume* FEEBoxLog   = caloBuildLogical(FEEBox,   FEEBoxMaterial, "caloFEEBoxLog",   0, G4Color::Yellow(),0,forceEdge);   
       G4LogicalVolume* FEEBoxInLog = caloBuildLogical(FEEBoxIn, vacuumMaterial, "caloFEEBoxInLog", 0, G4Color::Black(),0,0);   
       G4LogicalVolume* FEECardLog  = caloBuildLogical(FEECard,  FEEMaterial,    "caloFEECardROLog",isROVisible, G4Color::Green(),isROSolid,forceEdge);       

       // cards touch top of the box    
       double dFEEsize = FEEBoxHalfZ-0.5*FEEBoxThickness - FEEHalfZ;
       pv = new G4PVPlacement(0,G4ThreeVector(0, 0, -0.5*FEEBoxThickness), FEEBoxInLog,"caloFEEBoxIn_PV", FEEBoxLog, false,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(ROHalfX, 0, -dFEEsize),  FEECardLog,"caloFEECardPV_0", FEEBoxInLog, true,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(-ROHalfX, 0, -dFEEsize), FEECardLog,"caloFEECardPV_1", FEEBoxInLog, true,1,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

       // add sensitive detector    
       G4VSensitiveDetector* crCardSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloReadoutCard());
       if (crCardSD) FEECardLog->SetSensitiveDetector(crCardSD);


       //----------------------
       // Build backplate and readouts in the holes
       //----
       int nTotCrystal(0);
       for (int i=0;i<idisk;++i) nTotCrystal+=cal.disk(idisk).nCrystals();

       G4Tubs* backPlate = new G4Tubs("caloBackPlate",BPInnerRadius,BPOuterRadius,holeHalfZ,0,CLHEP::twopi);       
       G4LogicalVolume* backPlateLog = caloBuildLogical(backPlate, backPlateMaterial, "caloBackPlateLog", isROVisible, G4Color::Blue(),isROSolid,forceEdge);

       for(unsigned ic=0; ic <cal.disk(idisk).nCrystals(); ++ic)
       {	      
	     G4int id = nTotCrystal+ic;
	     std::ostringstream name;name<<"caloHolePV_" <<id;
             CLHEP::Hep3Vector unitPosition = cal.disk(idisk).crystal(ic).localPosition();
             unitPosition.setZ(0.0);

             pv = new G4PVPlacement(0,unitPosition,holeBackLog,name.str(),backPlateLog,false,id,false);
             doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       }

       
             
       // Approximate model used to avoid recording the length of each strip - good enough in our case
       // two kind of copper strips: cut by inner hole and full length
       // strips cut by inner hole have return pipes, and adjacent strips need to have the same length
       std::vector<double> stripY, stripX0, stripX1;
       for (int is=0;is<int(BPOuterRadius/2.0/wrapperHalfXY);is+=1)
       {
           double ycrystal = is*2.0*wrapperHalfXY;
           
           int idxIn(-1),idxOut(-1);
           for (unsigned i=0; i<stepsInY.size();++i)  if (abs(stepsInY[i]-ycrystal) < 0.1) idxIn=i;       
           for (unsigned i=0; i<stepsOutY.size();++i) if (abs(stepsOutY[i]-ycrystal) < 0.1) idxOut=i;
           
           double x0  = (idxIn > -1) ? stepsInX[idxIn] : 0.0;
           double x1  = (idxOut > -1) ? stepsOutX[idxOut] : 0.0;           
           if (idxIn > -1  && idxIn+1<nstepsInX)   x0 = std::min(x0,stepsInX[idxIn+1]);
           if (idxIn > -1  && idxIn+1==nstepsInX)  x0 = 0.0;
           if (idxOut > -1 && idxOut+1<nstepsOutX) x1 = std::max(x1,stepsOutX[idxOut+1]);
           
           if (x1-x0<1e-3) continue;
           stripY.push_back((2*is+1)*wrapperHalfXY);
           stripX0.push_back(x0);      
           stripX1.push_back(x1);      
       
           //hack to make it look more like the drawing for the current calorimeter         
           if (is<8 && is%2==1 )
           {
              double mi = std::min(stripX0[is],stripX0[is-1]);
              double ma = std::max(stripX1[is],stripX1[is-1]);
              stripX0[is]  = mi;
              stripX0[is-1]= mi;            
              stripX1[is]  = ma;
              stripX1[is-1]= ma;            
           }
       }
       
       //build the strips
       for (unsigned i=0;i<stripY.size();++i) 
       {                                
           double zpos = holeHalfZ-stripHalfZ;
           if (stripX0[i]>0.1)
           {
               G4Box* bstrip = new G4Box("strip",0.5*(stripX1[i]-stripX0[i]),stripHalfY, stripHalfZ);
               G4LogicalVolume* bstripLog = caloBuildLogical(bstrip, stripMaterial, "caloStripLog", isROVisible, G4Color::Red(),isROSolid,forceEdge);
               double xpos = stripX0[i]+0.5*(stripX1[i]-stripX0[i]);

               pv = new G4PVPlacement(0,G4ThreeVector(xpos,  stripY[i], zpos), bstripLog,"caloStripPV", backPlateLog, false,0,false);
               doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
  	       pv = new G4PVPlacement(0,G4ThreeVector(-xpos, stripY[i], zpos), bstripLog,"caloStripPV", backPlateLog, false,0,false);
               doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
               pv = new G4PVPlacement(0,G4ThreeVector(xpos, -stripY[i], zpos), bstripLog,"caloStripPV", backPlateLog, false,0,false);
               doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
  	       pv = new G4PVPlacement(0,G4ThreeVector(-xpos, -stripY[i], zpos),bstripLog,"caloStripPV", backPlateLog, false,0,false);
               doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
           } 
           else 
           {
               G4Box* bstrip = new G4Box("strip",stripX1[i],stripHalfY, stripHalfZ);
               G4LogicalVolume* bstripLog = caloBuildLogical(bstrip, stripMaterial, "caloStripLog", isROVisible, G4Color::Red(),isROSolid,forceEdge);

               pv = new G4PVPlacement(0,G4ThreeVector(0,  stripY[i], zpos), bstripLog,"caloStripPV", backPlateLog, false,0,false);
               doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
  	       pv = new G4PVPlacement(0,G4ThreeVector(0, -stripY[i], zpos), bstripLog,"caloStripPV", backPlateLog, false,0,false);
               doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
           }
       }


       //----------------------
       // Build FEE boxes and last cooling pipes
       //----                     
       G4Tubs* backPlateFEE = new G4Tubs("caloBackFEEPlate",BPInnerRadius,BPOuterRadius,BPFEEHalfZ,0,CLHEP::twopi);       
       G4LogicalVolume* backPlateFEELog = caloBuildLogical(backPlateFEE, vacuumMaterial, "caloBackPlateFEELog",0, G4Color::Green(),0,0);

       for(unsigned ic=0; ic <cal.disk(idisk).nCrystals(); ++ic)
       {	      
	   G4int id = nTotCrystal+ic;
	   std::ostringstream name;name<<"caloFEEPV_" <<id;
           CLHEP::Hep3Vector unitPosition = cal.disk(idisk).crystal(ic).localPosition();
           unitPosition.setZ(-BPFEEHalfZ + FEEBoxHalfZ); 

           pv = new G4PVPlacement(0,unitPosition,FEEBoxLog,name.str(),backPlateFEELog,false,id,false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       }

       //just add the cooling pipes outside, like that...
       //Trick: pipe in xy plane going from -alpha to beta = pipe going from pi-beta to pi+alpha, then rotated by pi around y axis
       G4Torus* BPPipe1 = new G4Torus("BPPipe1",BPPipeRadiusHigh-BPPipeThickness, BPPipeRadiusHigh,BPPipeTorRadiusHigh, 0.5*CLHEP::pi, 0.9*CLHEP::pi);
       G4Torus* BPPipe2 = new G4Torus("BPPipe2",BPPipeRadiusLow-BPPipeThickness,  BPPipeRadiusLow, BPPipeTorRadiusHigh, 0.5*CLHEP::pi, 0.7*CLHEP::pi);
       G4Torus* BPPipe3 = new G4Torus("BPPipe3",BPPipeRadiusHigh-BPPipeThickness, BPPipeRadiusHigh,BPPipeTorRadiusLow,  0.5*CLHEP::pi, 0.9*CLHEP::pi);
       G4Torus* BPPipe4 = new G4Torus("BPPipe4",BPPipeRadiusLow-BPPipeThickness,  BPPipeRadiusLow, BPPipeTorRadiusLow,  0.5*CLHEP::pi, 0.7*CLHEP::pi);

       G4LogicalVolume* BPPipe1Log = caloBuildLogical(BPPipe1, pipeMaterial, "caloBPPipe11Log",isROVisible,G4Color::Blue(),isROSolid,forceEdge);      
       G4LogicalVolume* BPPipe2Log = caloBuildLogical(BPPipe2, pipeMaterial, "caloBPPipe12Log",isROVisible,G4Color::Blue(),isROSolid,forceEdge);
       G4LogicalVolume* BPPipe3Log = caloBuildLogical(BPPipe3, pipeMaterial, "caloBPPipe13Log",isROVisible,G4Color::Red() ,isROSolid,forceEdge);      
       G4LogicalVolume* BPPipe4Log = caloBuildLogical(BPPipe4, pipeMaterial, "caloBPPipe14Log",isROVisible,G4Color::Red() ,isROSolid,forceEdge);

       G4RotationMatrix* rotY = new G4RotationMatrix(CLHEP::HepRotation::IDENTITY);rotY->rotateY(CLHEP::pi);

       pv = new G4PVPlacement(rotY,G4ThreeVector(0,0, BPFEEHalfZ- BPPipeRadiusHigh),BPPipe1Log,"caloBPPipe1PV",backPlateFEELog,false,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(0,0, BPFEEHalfZ- BPPipeRadiusHigh),BPPipe2Log,"caloBPPipe2PV",backPlateFEELog,false,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(0,0, BPFEEHalfZ- BPPipeRadiusHigh),BPPipe3Log,"caloBPPipe3PV",backPlateFEELog,false,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       pv = new G4PVPlacement(rotY,G4ThreeVector(0,0, BPFEEHalfZ- BPPipeRadiusHigh),BPPipe4Log,"caloBPPipe4PV",backPlateFEELog,false,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       


       //----------------------
       // Make the final full back plate
       //----                     
       G4Tubs* fullBackPlateFEE = new G4Tubs("caloFullBackPlate",BPInnerRadius,BPOuterRadius,BPFEEHalfZ+holeHalfZ,0,CLHEP::twopi);       
       G4LogicalVolume* fullBackPlateFEELog = caloBuildLogical(fullBackPlateFEE, vacuumMaterial, "caloFullBackPlateLog", 0, G4Color::Black(),0,0);

       pv = new G4PVPlacement(0,G4ThreeVector(0,0,-BPFEEHalfZ),backPlateLog,"caloFullBP1PV",fullBackPlateFEELog,false,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(0,0,holeHalfZ),backPlateFEELog,"caloFullBP2PV",fullBackPlateFEELog,false,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

       return fullBackPlateFEELog;
   }
  
  
  
  //--------------------------------------------------------------------------------------------------------------------------------
  // build crate
  G4LogicalVolume* caloBuildCrate(const SimpleConfig& config, MaterialFinder& materialFinder, const DiskCalorimeter& cal)
  {     
       const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
       geomOptions->loadEntry( config, "calorimeterCrate", "calorimeter.crate" );
       geomOptions->loadEntry( config, "calorimeterCrateBoard", "calorimeter.crateBoard" );

       const bool isCrateVisible       = geomOptions->isVisible("calorimeterCrate"); 
       const bool isCrateBoardVisible  = geomOptions->isVisible("calorimeterCrateBoard"); 
       const bool isCrateSolid         = geomOptions->isSolid("calorimeterCrate"); 
       const bool isCrateBoardSolid    = geomOptions->isSolid("calorimeterCrateBoard"); 
       const bool doSurfaceCheck       = config.getBool("g4.doSurfaceCheck",false) || config.getBool("calorimeter.doSurfaceCheck",false);
       const int  verbosityLevel       = config.getInt("calorimeter.verbosityLevel",1);
       const bool forceEdge            = config.getBool("g4.forceEdge",false);


       G4Material* vacuumMaterial        = materialFinder.get("calorimeter.vacuumMaterial");
       G4Material* crateMaterial         = materialFinder.get("calorimeter.crateMaterial");
       G4Material* crateBottomAMaterial  = materialFinder.get("calorimeter.crateMaterial");
       G4Material* crateShieldMaterial   = materialFinder.get("calorimeter.shieldMaterial");    
       G4Material* radiatorMaterial      = materialFinder.get("calorimeter.radiatorMaterial");
       G4Material* activeStripMaterial   = materialFinder.get("calorimeter.activeStripMaterial");
       G4Material* passiveStripMaterial  = materialFinder.get("calorimeter.passiveStripMaterial");
 

       G4double crystalHalfZ             = cal.caloInfo().getDouble("crystalZLength")/2.0;    
       G4int nBoards                     = cal.caloInfo().getInt("numberOfBoards");
       G4double crateDX                  = cal.caloInfo().getDouble("crateXLength")/2.0;
       G4double crateDY                  = cal.caloInfo().getDouble("crateYLength")/2.0;
       G4double crateDZ                  = cal.caloInfo().getDouble("crateZLength")/2.0;
       G4double crateADZ                 = crystalHalfZ;
       G4double crateFShieldDisp         = cal.caloInfo().getDouble("crateFShieldDeltaZ");
       G4double crateFShieldThick        = cal.caloInfo().getDouble("crateFShieldThickness");
       G4double crateBottomThick         = cal.caloInfo().getDouble("crateBShieldThickness");
       G4double crateTopThick            = cal.caloInfo().getDouble("crateTThickness");
       G4double crateSideThick           = cal.caloInfo().getDouble("crateSThickness");
       G4double crateFShieldDY           = cal.caloInfo().getDouble("crateFShieldYLength")/2.0;
       G4double crateFullDZ              = crateDZ+crateFShieldDisp/2.0+crateFShieldThick/2.0;
       G4double crateFullDY              = crateDY+crateBottomThick/2.0;
       G4double crateBoxInDY             = crateDY-crateTopThick/2.0;
       G4double crateBDZ                 = crateFullDZ - crateADZ;

       G4double radiatorDY               = cal.caloInfo().getDouble("radiatorThickness")/2.0;
       G4double activeStripDY            = cal.caloInfo().getDouble("activeStripThickness")/2.0;
       G4double passiveStripDY           = cal.caloInfo().getDouble("passiveStripThickness")/2.0;


       G4VPhysicalVolume* pv;
       
       // define the crate box and shield pieces
       G4Box *crateFullBox = new G4Box("ccrateFullBox",crateDX,crateFullDY,crateFullDZ);
       G4Box *crateBox     = new G4Box("ccrateBoxTop", crateDX,crateDY,crateDZ);
       G4Box *crateBoxIn   = new G4Box("ccrateBoxSide",crateDX-crateSideThick,crateBoxInDY,crateDZ);
       G4Box *crateBottomA = new G4Box("ccrateBottomA",crateDX,crateBottomThick/2.0,crateADZ);
       G4Box *crateBottomB = new G4Box("ccrateBottomB",crateDX,crateBottomThick/2.0,crateBDZ);
       G4Box *crateShieldF = new G4Box("ccrateShieldF",crateDX,crateFShieldDY,crateFShieldThick/2.0);
       
       //define the logical volumes
       G4LogicalVolume *crateFullBoxLog = caloBuildLogical(crateFullBox, vacuumMaterial,      "ccrateFullBoxLog",0,G4Color::White(),0,0);   
       G4LogicalVolume *crateBoxLog     = caloBuildLogical(crateBox,     crateMaterial,       "ccrateBoxLog",isCrateVisible,G4Color::Blue(),isCrateSolid,forceEdge);   
       G4LogicalVolume *crateBoxInLog   = caloBuildLogical(crateBoxIn,   vacuumMaterial,      "ccrateBoxInLog",isCrateVisible,G4Color::Green(),isCrateSolid,forceEdge);   
       G4LogicalVolume *crateBottomALog = caloBuildLogical(crateBottomA, crateBottomAMaterial,"ccrateBottomALog",isCrateVisible,G4Color::Red(),isCrateSolid,forceEdge);   
       G4LogicalVolume *crateBottomBLog = caloBuildLogical(crateBottomB, crateMaterial,       "ccrateBottomBLog",isCrateVisible,G4Color::Magenta(),isCrateSolid,forceEdge);   
       G4LogicalVolume *crateShieldFLog = caloBuildLogical(crateShieldF, crateShieldMaterial, "ccrateShieldFLog",isCrateVisible,G4Color::Yellow(),isCrateSolid,forceEdge);   
              
       //place the crate box and shield pieces
       G4double crateFZpos  = -crateFullDZ + crateFShieldThick/2.0;
       G4double crateAZpos  = -crateFullDZ + crateADZ;
       G4double crateBZpos  = -crateFullDZ + 2.0*crateADZ + crateBDZ;
       G4double crateABYpos = -crateFullDY + crateBottomThick/2.0;
       G4double crateFYpos  = -crateFullDY + crateBottomThick + crateFShieldDY;

       pv = new G4PVPlacement(0,G4ThreeVector(0.0,crateBottomThick/2.0,crateFullDZ-crateDZ),crateBoxLog,"ccrateBoxPV",crateFullBoxLog,false,0,false);
       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,-crateTopThick/2.0,0),crateBoxInLog,"ccrateBoxInPV",crateBoxLog,false,0,false);
       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);       
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,crateABYpos,crateAZpos),crateBottomALog,"ccrateBottomAPV",crateFullBoxLog,false,0,false);
       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,crateABYpos,crateBZpos),crateBottomBLog,"ccrateBottomBPV",crateFullBoxLog,false,0,false);
       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
       if (crateFShieldThick>0.01) pv = new G4PVPlacement(0,G4ThreeVector(0.0, crateFYpos,crateFZpos),crateShieldFLog,"ccrateShieldFPV",crateFullBoxLog,false,0,false);
       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);

      
       // add the boards
       G4double boardDX          = crateDX-crateSideThick;    
       G4double boardDY          = radiatorDY+activeStripDY+passiveStripDY;
       G4double boardDZ          = crateDZ;
       G4double radiatorPosY     = -boardDY + radiatorDY;   
       G4double activeStripPosY  = radiatorPosY+radiatorDY+activeStripDY;
       G4double passiveStripPosY = activeStripPosY + activeStripDY + passiveStripDY;
       
       G4Box* boardCrate        = new G4Box("ccrateBoard",boardDX,boardDY,boardDZ);
       G4Box* radiatorBoard     = new G4Box("ccrateRadiator",boardDX,radiatorDY,boardDZ);
       G4Box* activeStripBoard  = new G4Box("ccrateActiveStrip",boardDX,activeStripDY,boardDZ);
       G4Box* passiveStripBoard = new G4Box("ccratePassiveStrip",boardDX,passiveStripDY,boardDZ);
      
       G4LogicalVolume* boardCrateLog        = caloBuildLogical(boardCrate, vacuumMaterial, "ccrateBoardLog",isCrateBoardVisible,G4Color::Green(),0,0);   
       G4LogicalVolume* radiatorBoardLog     = caloBuildLogical(radiatorBoard,radiatorMaterial, "ccrateRadiatorLog",0,G4Color::Black(),isCrateBoardSolid,forceEdge);   
       G4LogicalVolume* activeStripBoardLog  = caloBuildLogical(activeStripBoard, activeStripMaterial, "ccrateActiveStripLog",0,G4Color::Black(),isCrateBoardSolid,forceEdge);   
       G4LogicalVolume* passiveStripBoardLog = caloBuildLogical(passiveStripBoard, passiveStripMaterial, "ccratePassiveStripLog",0,G4Color::Black(),isCrateBoardSolid,forceEdge);   

       pv = new G4PVPlacement(0,G4ThreeVector(0.0,radiatorPosY,0.0),radiatorBoardLog,"ccrateRadiatorPV",boardCrateLog,false,0,false);
       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,activeStripPosY,0.0),activeStripBoardLog,"ccrateActiveStripPV",boardCrateLog,false,0,false);
       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,passiveStripPosY,0.0),passiveStripBoardLog,"ccratePassiveStripPV",boardCrateLog,false,0,false);
       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);

       // put boards in a single crate
       G4double deltaY  = 2.0*crateBoxInDY/float(nBoards+1);

       for (G4int ibrd=0; ibrd < nBoards; ++ibrd)
       {
	   std::ostringstream boardPV; boardPV<<"ccrateBoardPV_" <<ibrd;
	   pv = new G4PVPlacement(0,G4ThreeVector(0.0,(ibrd+1)*deltaY-crateBoxInDY,0), boardCrateLog, boardPV.str(), crateBoxInLog, false, ibrd,false);
	   doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
	   
           G4VSensitiveDetector* crCrate = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloCrate());
	   if (crCrate) activeStripBoardLog->SetSensitiveDetector(crCrate);	   
       }
      
      return crateFullBoxLog;
   }

  
  //--------------------------------------------------------------------------------------------------------------------------------
  // build full FEB
  G4LogicalVolume* caloBuildFEB(const SimpleConfig& config, MaterialFinder& materialFinder, const DiskCalorimeter& cal)
  {                  
       const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
       geomOptions->loadEntry( config, "calorimeterCrate", "calorimeter.crate" );

       const bool isCrateVisible  = geomOptions->isVisible("calorimeterCrate"); 
       //const bool isCrateVisible  = config.getBool("calorimeter.crateVisible",false);
       const bool doSurfaceCheck  = config.getBool("g4.doSurfaceCheck",false) || config.getBool("calorimeter.doSurfaceCheck",false);
       const int  verbosityLevel  = config.getInt("calorimeter.verbosityLevel",1);

       G4Material* vacuumMaterial = materialFinder.get("calorimeter.vacuumMaterial");
       G4VPhysicalVolume* pv;
       
       
       // get me a lil' crate
       G4LogicalVolume* crateBoxLog = caloBuildCrate(config, materialFinder, cal);
       G4Box* crate = static_cast<G4Box*>(crateBoxLog->GetSolid());
 
       G4int nCrates            = cal.caloInfo().getInt("numberOfCrates");
       G4int nCratesBeforeSpace = cal.caloInfo().getInt("nCrateBeforeSpace");
       G4double phi0Crate       = cal.caloInfo().getDouble("cratephi0")*CLHEP::degree;
       G4double deltaPhiCrate   = cal.caloInfo().getDouble("crateDeltaPhi")*CLHEP::degree; 
       G4double crateRadIn      = cal.caloInfo().getDouble("crateInnerRadius");
       G4double crateRadOut     = crateRadIn+1.005*sqrt(4.0*crate->GetYHalfLength()*crate->GetYHalfLength()+crate->GetXHalfLength()*crate->GetXHalfLength());
       G4double crateHalfLength = crate->GetZHalfLength();   
       G4double cratePosY       = crateRadIn+1.001*crate->GetYHalfLength();    
       
      
       G4Tubs* calorimeterFEB = new G4Tubs("caloFEB",crateRadIn, crateRadOut,crateHalfLength,-phi0Crate, CLHEP::pi+2*phi0Crate);       
       G4LogicalVolume* calorimeterFEBLog  = caloBuildLogical(calorimeterFEB, vacuumMaterial, "caloFEBLog",0,G4Color::Black(),0,0);   
  
       G4RotationMatrix rotCrate = G4RotationMatrix();
       G4double phiCrate(0);       
       for (G4int icrt=0;icrt < nCrates;++icrt)
       {		   
           if (icrt<nCrates/2) 
           {
              phiCrate = CLHEP::pi/2. - deltaPhiCrate*(icrt+1);
              if (icrt >= nCratesBeforeSpace) phiCrate -= deltaPhiCrate;
           }
           else 
	   {
	       phiCrate = CLHEP::pi/2. + deltaPhiCrate*(icrt-nCrates/2+1);
	       if (icrt >= nCrates/2+nCratesBeforeSpace) phiCrate += deltaPhiCrate;		
	   }
           
           if (icrt==nCrates/2) rotCrate.rotateZ(CLHEP::pi/2.+deltaPhiCrate+deltaPhiCrate/3.0);
	   else if (icrt<nCrates/2) rotCrate.rotateZ(-deltaPhiCrate);
	   else rotCrate.rotateZ(deltaPhiCrate);                     
           if (icrt == nCratesBeforeSpace) rotCrate.rotateZ(-deltaPhiCrate);
           if (icrt == nCrates/2+nCratesBeforeSpace) rotCrate.rotateZ(deltaPhiCrate);
           	   
           G4ThreeVector posCrate   = G4ThreeVector(cratePosY*std::cos(phiCrate),cratePosY*std::sin(phiCrate),0.0);
	   G4Transform3D crateCoord = G4Transform3D(rotCrate,posCrate);

	   std::ostringstream cratePV; cratePV<<"caloCratePV_" <<icrt;
	   pv = new G4PVPlacement(crateCoord, crateBoxLog, cratePV.str(), calorimeterFEBLog, false, icrt, false);
	   doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);	
        }
        
 	if (config.getBool("ds.hasCableRunCal",false))
        {
             double crRin  = config.getDouble("ds.CableRunCal.Rin")*CLHEP::mm;
	     double crRout = config.getDouble("ds.CableRunCal.Rout")*CLHEP::mm;
	     double phi0   = config.getDouble("ds.CableRunCal.phi0")*CLHEP::degree;
	     double dPhi   = config.getDouble("ds.CableRunCal.dPhi")*CLHEP::degree;

             G4Material* cableMaterial = findMaterialOrThrow(config.getString("ds.CableRunCal.material"));

	     G4Tubs* ccrTub = new G4Tubs("caloCableRunCalTub",crRin, crRout,crateHalfLength - 5.0,phi0,dPhi);
             G4LogicalVolume* ccrTubLog = caloBuildLogical(ccrTub, cableMaterial, "caloCableRunCalTubLog",isCrateVisible,G4Color::Yellow(),0,0);

	     G4RotationMatrix ccrRot = G4RotationMatrix();
	     CLHEP::Hep3Vector calCableRunLoc(0.0,0.0,0.0);
	     G4Transform3D ccrCoord = G4Transform3D(ccrRot,calCableRunLoc);

	     pv = new G4PVPlacement(ccrCoord,ccrTubLog,"caloCableRunCalTub_PV",calorimeterFEBLog, false, 0, false);
	     doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

         }         
 	 if ( config.getBool("ds.hasCableRunTrk",false))
         {
	     double crRin = config.getDouble("ds.CableRunTrk.Rin")*CLHEP::mm;
	     double crRout = config.getDouble("ds.CableRunTrk.Rout")*CLHEP::mm;
	     double phi01 = config.getDouble("ds.CableRunTrk.phi0")*CLHEP::degree;
	     double dPhi = config.getDouble("ds.CableRunTrk.dPhi")*CLHEP::degree;
	     double phi02 = 180.0*CLHEP::degree - phi01 - dPhi;
	     if ( phi01 + dPhi > 180.0*CLHEP::degree + phi0Crate )
             {
		 dPhi = 179.5*CLHEP::degree + phi0Crate - phi01;
		 phi02 = 180.0*CLHEP::degree - phi01 - dPhi;
	     }

             G4Material* cableMaterial = findMaterialOrThrow(config.getString("ds.CableRunTrk.material"));

	     G4Tubs* ccr1Tub = new G4Tubs("caloCableRun1TrkTub",crRin, crRout, crateHalfLength - 5.0,phi01, dPhi);
	     G4Tubs* ccr2Tub = new G4Tubs("caloCableRun2TrkTub",crRin, crRout,crateHalfLength - 5.0,phi02, dPhi);
             G4LogicalVolume* ccr1TubLog = caloBuildLogical(ccr1Tub, cableMaterial, "caloccr1TubLog",isCrateVisible,G4Color::Yellow(),0,0);
             G4LogicalVolume* ccr2TubLog = caloBuildLogical(ccr2Tub, cableMaterial, "caloccr2TubLog",isCrateVisible,G4Color::Yellow(),0,0);

	     CLHEP::Hep3Vector trkCableRunLoc(0.0,0.0,0.0);
	     G4RotationMatrix ccrRot = G4RotationMatrix();
	     G4Transform3D ccrCoord = G4Transform3D(ccrRot,trkCableRunLoc);

	     pv = new G4PVPlacement(ccrCoord,ccr1TubLog,"TrkCableRun1InCalFeb",calorimeterFEBLog,false, 0, false);
	     doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

	     pv = new G4PVPlacement(ccrCoord,ccr2TubLog,"TrkCableRun2InCalFeb",calorimeterFEBLog,false, 0, false);
	     doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);			     
	 }
         
         return calorimeterFEBLog;
     }




     //--------------------------------------------------------------------------------------------------------------------------------
     // utility for Logical volume
     G4LogicalVolume* caloBuildLogical(G4VSolid* solid, G4Material* mat, const G4String& name, bool isVisible, const G4Color& color, bool isSolid, bool forceEdge)
     {
        G4LogicalVolume* logical = new G4LogicalVolume(solid, mat, name);
        G4VisAttributes* visAtt = new G4VisAttributes(isVisible, color);
        if (!isVisible) logical->SetVisAttributes(G4VisAttributes::Invisible);
        else 
        {
          visAtt->SetForceSolid(isSolid);
          visAtt->SetForceAuxEdgeVisible(forceEdge);
          logical->SetVisAttributes(visAtt);
        }
        return logical;        
     } 



}
