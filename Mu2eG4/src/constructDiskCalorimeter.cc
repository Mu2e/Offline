//
// Free function to create the Disk calorimeter - please keep it clean
//
//
// Original author Bertrand Echenard
//
//
// Note about chimes in the disk. There are two ways to add the chimes between the crystals and the disk edges. 
//      1) create each edge separately (by subtracting the inner hole or as the intersection between the outside edge and the bar)
//      2) make an union comprising of all the inner chimes described as boxes, then subtract the hole inside. For the outer chimes
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

       G4PVPlacement *pv;


    //--------------------------------------
    // Construct calorimeter mother volume

       const DiskCalorimeter& cal  = *(GeomHandle<DiskCalorimeter>());
       
       const int  crateVersion     = config.getInt("calorimeter.crateVersion",1);
       const unsigned int nDisks   = cal.nDisk();
       G4double mother_inRadius    = cal.caloInfo().envelopeRadiusIn();
       G4double mother_outRadius   = cal.caloInfo().envelopeRadiusOut();
       G4double mother_z0          = cal.caloInfo().envelopeZ0();
       G4double mother_z1          = cal.caloInfo().envelopeZ1();
       G4double mother_zlength     = mother_z1-mother_z0;
       G4double mother_zCenter     = (mother_z1+mother_z0)/2.0;

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
             G4Tubs* backPlate  = (backPlateLog != nullptr) ? static_cast<G4Tubs*>(backPlateLog->GetSolid()) : nullptr; 
             G4Tubs* disk       = static_cast<G4Tubs*>(diskLog->GetSolid());

             G4double R0disk    = disk->GetInnerRadius();
             G4double R1disk    = disk->GetOuterRadius();
             G4double zHalfDisk = disk->GetZHalfLength();
             G4double zHalfFP   = (frontPlateLog != nullptr) ? frontPlate->GetZHalfLength() : 0;
             G4double zHalfBP   = (backPlateLog != nullptr) ? backPlate->GetZHalfLength() : 0;
             G4double zHalftot  = zHalfFP+zHalfDisk+zHalfBP;
             G4double vdThick(1.0*CLHEP::mm);
            
             double diskpar[5] = {R0disk-vdThick,R1disk,zHalftot+0.5*vdThick,0,CLHEP::twopi};
             std::ostringstream discname;  discname<<"caloDisk_" <<idisk;

	     //origin gives the position of the center of the disk, irrespective of the coordinate origin set in the calo description
	     G4ThreeVector posDisk = cal.disk(idisk).geomInfo().origin() - posCaloMother;

	     calorimeterDisk[idisk] =  nestTubs(discname.str(),diskpar,vacuumMaterial,&cal.disk(idisk).geomInfo().rotation(),posDisk,
                                                calorimeterInfo,idisk,
                        		        isCalorimeterVisible,G4Colour::White(),0,forceEdge,true,doSurfaceCheck );
                                                
             if (frontPlateLog!=nullptr) pv = new G4PVPlacement(0,G4ThreeVector(0,0,-zHalftot+zHalfFP),frontPlateLog,"caloFP_PV",calorimeterDisk[idisk].logical,false,0,false);
             if (frontPlateLog!=nullptr) doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
             
             pv = new G4PVPlacement(0,G4ThreeVector(0,0,-zHalftot+2*zHalfFP+zHalfDisk),diskLog,"caloCase_PV",calorimeterDisk[idisk].logical,false,0,false);
             doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
             
             if (backPlateLog!=nullptr) pv = new G4PVPlacement(0,G4ThreeVector(0,0,+zHalftot-zHalfBP),backPlateLog,"caloBP_PV",calorimeterDisk[idisk].logical,false,0,false);
             if (backPlateLog!=nullptr) doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
             


	     if ( crateVersion > 1 )
             {
   	        std::ostringstream cratename; cratename<<"caloFEB_" <<idisk;
                double FEBpar[5] = {FEB->GetInnerRadius()-vdThick,FEB->GetOuterRadius()+vdThick,FEB->GetZHalfLength()+0.5*vdThick,FEB->GetStartPhiAngle(),FEB->GetDeltaPhiAngle()};
                                                
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
  // Front plate: a styrofoam panel with the pipes inside, surrounded by two cooling pipes.
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
       G4Material* FPStyrofoamMaterial   = materialFinder.get("calorimeter.FPStyrofoamMaterial");
       G4Material* FPCarbonMaterial      = materialFinder.get("calorimeter.FPCarbonMaterial");
       G4Material* pipeMaterial          = materialFinder.get("calorimeter.pipeMaterial");

       G4double diskCaseOuterRadius      = cal.geomInfo().diskCaseRadiusOut().at(idisk);
       G4double outerCaseThickness       = cal.caloInfo().caseThicknessOut();
       G4double outerRingRadius          = diskCaseOuterRadius + outerCaseThickness;
       G4double FPInnerRadius            = cal.caloInfo().FPInnerRadius();
       G4double FPStyroHalfThick         = cal.caloInfo().FPStyrofoamZLength()/2.0;  
       G4double FPCarbonHalfThick        = cal.caloInfo().FPCarbonZLength()/2.0;  
       G4double frontPanelHalfThick      = FPStyroHalfThick+2.0*FPCarbonHalfThick;

       G4int nPipes                      = cal.caloInfo().nPipes();      
       G4double pipeRadius               = cal.caloInfo().pipeRadius();
       G4double pipeThickness            = cal.caloInfo().pipeThickness();
       std::vector<double> pipeTorRadius = cal.caloInfo().pipeTorRadius();
       G4double pipeSeparation           = cal.caloInfo().pipeSeparation();    
       G4double coolingPipeRadius        = cal.caloInfo().coolFPPipeRadius();
       G4double pipeZpos                 = -FPStyroHalfThick+pipeRadius;
       G4double coolpipeZpos             = -frontPanelHalfThick+2.0*FPCarbonHalfThick+pipeRadius;
       

       if (nPipes==0) return nullptr;

       //this is the full front panel
       G4Tubs* frontPlate = new G4Tubs("caloFrontPlate",FPInnerRadius,outerRingRadius+2*coolingPipeRadius,frontPanelHalfThick,0,CLHEP::twopi);       
       G4LogicalVolume* frontPlateLog = caloBuildLogical(frontPlate, vacuumMaterial, "caloFrontPlateLog",0,G4Color::White(),0,0);

       //styrofoam panel
       G4Tubs* frontPanelStyro = new G4Tubs("caloFPStyro",FPInnerRadius,outerRingRadius,FPStyroHalfThick,0,CLHEP::twopi);       
       G4LogicalVolume* frontPanelStyroLog = caloBuildLogical(frontPanelStyro, FPStyrofoamMaterial, "caloFPStyroLog",isPipeVisible,G4Color::Brown(),0,forceEdge);
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,0.0), frontPanelStyroLog, "caloFPStyroPV", frontPlateLog, false, 0, false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);                

       //carbon fiber panel
       G4Tubs* frontPanelCarb = new G4Tubs("caloFPCarb",FPInnerRadius,outerRingRadius,FPCarbonHalfThick,0,CLHEP::twopi);       
       G4LogicalVolume* frontPanelCarbLog = caloBuildLogical(frontPanelCarb, FPCarbonMaterial, "caloFPCarbLog",isPipeVisible,G4Color::Grey(),0,forceEdge);
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,-frontPanelHalfThick+FPCarbonHalfThick), frontPanelCarbLog, "caloFPCarbPV1", frontPlateLog, false, 0, false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);                
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,frontPanelHalfThick-FPCarbonHalfThick), frontPanelCarbLog, "caloFPCarbPV2", frontPlateLog, false, 0, false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);                
       
       //cooling pipes on the edge
       double angMax = std::acos((nPipes-0.5)*pipeSeparation/outerRingRadius)-0.1;
       G4Torus* coolFR = new G4Torus("caloCoolFP",coolingPipeRadius-pipeThickness, coolingPipeRadius, outerRingRadius+coolingPipeRadius, 0, CLHEP::twopi-angMax);
       G4LogicalVolume* coolFRLog  = caloBuildLogical(coolFR, pipeMaterial, "caloCoolFPLog",isPipeVisible,G4Color::Red(),isPipeSolid,0);
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,coolpipeZpos), coolFRLog, "caloCoolFPPV", frontPlateLog, false, 0, false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);                
       

       G4RotationMatrix* rotPipe1 = new G4RotationMatrix(CLHEP::HepRotation::IDENTITY);
       rotPipe1->rotateZ(CLHEP::pi/2.0);
       G4RotationMatrix* rotPipe2 = new G4RotationMatrix(CLHEP::HepRotation::IDENTITY);
       rotPipe2->rotateZ(1.5*CLHEP::pi);
       G4RotationMatrix* rotPipeFlat = new G4RotationMatrix(CLHEP::HepRotation::IDENTITY);
       rotPipeFlat->rotateX(CLHEP::pi/2.0);

       for (int ipipe=0; ipipe<nPipes; ++ipipe)
       {
           double xpipe  = (0.5+ipipe)*pipeSeparation;
           double angle  = std::asin(xpipe/pipeTorRadius[ipipe]); //angle taken w.r.t y axis!
           double length = sqrt(outerRingRadius*outerRingRadius-xpipe*xpipe) - (pipeTorRadius[ipipe]+pipeRadius)*cos(angle) - 2.0*pipeRadius;
           double y0     = (pipeTorRadius[ipipe]+pipeRadius)*cos(angle)+0.5*length;
           double y1     = -(pipeTorRadius[ipipe]+pipeRadius)*cos(angle)-0.5*length;

           G4Torus* pipe1 = new G4Torus("caloPipe",pipeRadius-pipeThickness, pipeRadius, pipeTorRadius[ipipe],angle, CLHEP::pi-2*angle);
           G4Tubs*  pipe2 = new G4Tubs("caloPipe2",pipeRadius-pipeThickness, pipeRadius,0.5*length,0,CLHEP::twopi);       
           G4LogicalVolume* pipe1Log = caloBuildLogical(pipe1, pipeMaterial, "caloPipe1Log",isPipeVisible,G4Color::Cyan(),isPipeSolid,forceEdge);
           G4LogicalVolume* pipe2Log = caloBuildLogical(pipe2, pipeMaterial, "caloPipe2Log",isPipeVisible,G4Color::Cyan(),isPipeSolid,forceEdge);

           pv = new G4PVPlacement(rotPipe1,G4ThreeVector(0,0,pipeZpos), pipe1Log, "caloPipePV", frontPanelStyroLog, false, ipipe, false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);                
           pv = new G4PVPlacement(rotPipe2,G4ThreeVector(0,0,pipeZpos), pipe1Log, "caloPipePV", frontPanelStyroLog, false, ipipe, false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);                                
           pv = new G4PVPlacement(rotPipeFlat,G4ThreeVector(xpipe, y0,pipeZpos), pipe2Log, "caloPipePV", frontPanelStyroLog, false, ipipe , false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);           
           pv = new G4PVPlacement(rotPipeFlat,G4ThreeVector(xpipe,-y0,pipeZpos), pipe2Log, "caloPipePV", frontPanelStyroLog, false, ipipe , false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);           
           pv = new G4PVPlacement(rotPipeFlat,G4ThreeVector(-xpipe,y1,pipeZpos), pipe2Log, "caloPipePV", frontPanelStyroLog, false, ipipe , false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);           
           pv = new G4PVPlacement(rotPipeFlat,G4ThreeVector(-xpipe,-y1,pipeZpos), pipe2Log, "caloPipePV", frontPanelStyroLog, false, ipipe , false);
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

       const bool isDiskVisible        = geomOptions->isVisible("calorimeterCase"); 
       const bool isDiskSolid          = geomOptions->isSolid("calorimeterCase"); 
       const bool isCrystalVisible     = geomOptions->isVisible("calorimeterCrystal"); 
       const bool isCrystalSolid       = geomOptions->isSolid("calorimeterCrystal");        
       const bool forceEdge            = config.getBool("g4.forceEdge",false);
       const bool doSurfaceCheck       = config.getBool("g4.doSurfaceCheck",false) || config.getBool("calorimeter.doSurfaceCheck",false);
       const int  verbosityLevel       = config.getInt("calorimeter.verbosityLevel",1);

       G4Material* vacuumMaterial      = materialFinder.get("calorimeter.vacuumMaterial");
       G4Material* shimeMaterial       = materialFinder.get("calorimeter.shimeMaterial");
       G4Material* crysMaterial        = materialFinder.get("calorimeter.crystalMaterial");
       G4Material* wrapMaterial        = materialFinder.get("calorimeter.wrapperMaterial");    
       G4Material* innerRingMaterial   = materialFinder.get("calorimeter.innerRingMaterial");
       G4Material* outerRingMaterial   = materialFinder.get("calorimeter.outerRingMaterial");


       G4double crystalHalfXY          = cal.caloInfo().crystalXYLength()/2.0;
       G4double crystalHalfZ           = cal.caloInfo().crystalZLength()/2.0;    
       G4double wrapperHalfThick       = cal.caloInfo().wrapperThickness()/2.0;    
       G4double wrapperHalfXY          = crystalHalfXY + 2*wrapperHalfThick;
       G4double wrapperHalfZ           = crystalHalfZ  + wrapperHalfThick;

       G4double diskCaseInnerRadius    = cal.geomInfo().diskCaseRadiusIn().at(idisk);
       G4double diskCaseOuterRadius    = cal.geomInfo().diskCaseRadiusOut().at(idisk);
       G4double innerCaseThickness     = cal.caloInfo().caseThicknessIn();
       G4double outerCaseThickness     = cal.caloInfo().caseThicknessOut();
       G4double outerRingEdgeDepth     = cal.caloInfo().outerRingEdgeZLength();
       G4double outerRingEdgeThickness = cal.caloInfo().outerRingEdgeRLength();     
       G4double outerCaseRadius        = diskCaseOuterRadius + outerCaseThickness;
       G4double innerCaseRadius        = diskCaseInnerRadius - innerCaseThickness;
       G4double diskDepth              = wrapperHalfZ;

       std::vector<double> chimesInX   = cal.caloInfo().chimesInsideX();
       std::vector<double> chimesInY   = cal.caloInfo().chimesInsideY();
       std::vector<double> chimesOutX  = cal.caloInfo().chimesOutsideX();
       std::vector<double> chimesOutY  = cal.caloInfo().chimesOutsideY();

       G4VPhysicalVolume* pv;


       //-----------------------------------------
       // Build wrapper crystal unit
       //----
       G4Box* crystal     = new G4Box("caloCrystal",    crystalHalfXY,crystalHalfXY,crystalHalfZ);
       G4Box* crystalWrap = new G4Box("caloCrystalWrap",wrapperHalfXY,wrapperHalfXY,wrapperHalfZ);

       G4LogicalVolume *crystalLog = caloBuildLogical(crystal,     crysMaterial, "caloCrystalLog",0, G4Color::Cyan(),0,0);
       G4LogicalVolume *wrapperLog = caloBuildLogical(crystalWrap, wrapMaterial, "caloCrystalWrapLog",isCrystalVisible,G4Color::Cyan(),isCrystalSolid,forceEdge);  

       // crystal center in z = half wrapper size w.r.t wrapper center since there is no wrapper at the back of the crystal
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,wrapperHalfThick),crystalLog,"caloCrysPV",wrapperLog,false,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

       G4VSensitiveDetector* ccSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloCrystal());
       if (ccSD) crystalLog->SetSensitiveDetector(ccSD);


       //------------------------------------------------------------
       // Build disk cases (inner, crystal and outer cases + two external rings)
       //----
       G4Tubs* fullCaseDisk    = new G4Tubs("caloDiskFullCase",   innerCaseRadius,outerCaseRadius+outerRingEdgeThickness,diskDepth,0,CLHEP::twopi);       
       G4Tubs* crystalCaseDisk = new G4Tubs("caloDiskCrystalCase",innerCaseRadius,outerCaseRadius,                       diskDepth,0,CLHEP::twopi);       
       G4Tubs* ringCaseDisk    = new G4Tubs("caloDiskRingCase",   outerCaseRadius,outerCaseRadius+outerRingEdgeThickness,outerRingEdgeDepth,0,CLHEP::twopi);       

       G4LogicalVolume* fullCaseDiskLog    = caloBuildLogical(fullCaseDisk,    vacuumMaterial,    "caloDiskFullCaseLog",   0, G4Color::Black(),0,0);   
       G4LogicalVolume* crystalCaseDiskLog = caloBuildLogical(crystalCaseDisk, shimeMaterial,     "caloDiskCrystalCaseLog",isDiskVisible,G4Color::Yellow(),isDiskSolid,forceEdge);   
       G4LogicalVolume* ringCaseDiskLog    = caloBuildLogical(ringCaseDisk,    outerRingMaterial, "caloDiskRingCaseLog",   isDiskVisible,G4Color::Yellow(),isDiskSolid,forceEdge);   

       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0,0.0),crystalCaseDiskLog,"caloDiskCrystalCasePV",fullCaseDiskLog,true,0,false);
       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0,diskDepth-outerRingEdgeDepth),ringCaseDiskLog,"caloDiskRingCasePV",fullCaseDiskLog,true,0,false);
       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(0.0,0,-diskDepth+outerRingEdgeDepth),ringCaseDiskLog,"caloDiskRingCasePV",fullCaseDiskLog,true,0,false);
       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);




       //--------------------------------------------------------------------
       // add the chimes in the inner / outer part of the disk (see top Note). This is surprisingly compact...
       G4Tubs* holeIn = new G4Tubs("InnerDiskEdge",0,innerCaseRadius,diskDepth,0,CLHEP::twopi);           
       for (unsigned i=0;i<chimesInX.size();++i)
       {            
           G4Box* box = new G4Box("ChimeIn",chimesInX[i],wrapperHalfXY,crystalHalfZ);
           G4SubtractionSolid* cutChime = new G4SubtractionSolid("Cropped chime", box, holeIn,0, G4ThreeVector(0,-chimesInY[i],0));
           G4LogicalVolume* cutChimeLog = caloBuildLogical(cutChime, innerRingMaterial, "cutChimeLog",  isDiskVisible, G4Color::Yellow(), isCrystalSolid, forceEdge);   
           pv = new G4PVPlacement(0,CLHEP::Hep3Vector(0,chimesInY[i],0),cutChimeLog,"Chimes",crystalCaseDiskLog,true,0,false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       }

       for (unsigned i=0;i<chimesOutX.size();++i)
       {
           double xlen = sqrt(outerCaseRadius*outerCaseRadius-chimesOutY[i]*chimesOutY[i])-chimesOutX[i];
           G4Box* box = new G4Box("ChimeOut",xlen,wrapperHalfXY,crystalHalfZ);
           G4IntersectionSolid* cutChime1 = new G4IntersectionSolid("Cropped chime", crystalCaseDisk, box,0, G4ThreeVector(chimesOutX[i]+xlen,chimesOutY[i],0));
           G4LogicalVolume* cutChime1Log = caloBuildLogical(cutChime1, outerRingMaterial, "cutChimeLog",  isDiskVisible, G4Color::Yellow(), isCrystalSolid, forceEdge);   
           pv = new G4PVPlacement(0,CLHEP::Hep3Vector(0,0,0),cutChime1Log,"Chimes",crystalCaseDiskLog,true,0,false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

           G4IntersectionSolid* cutChime2 = new G4IntersectionSolid("Cropped chime", crystalCaseDisk, box,0, G4ThreeVector(-chimesOutX[i]-xlen,chimesOutY[i],0));
           G4LogicalVolume* cutChime2Log = caloBuildLogical(cutChime2, outerRingMaterial, "cutChimeLog",  isDiskVisible, G4Color::Yellow(), isCrystalSolid, forceEdge);   
           pv = new G4PVPlacement(0,CLHEP::Hep3Vector(0,0,0),cutChime2Log,"Chimes",crystalCaseDiskLog,true,0,false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       }

       //add the two pieces at the top and bottom
       if (chimesOutY.size()>0)
       {
           double ylow = *std::min_element(chimesOutY.begin(),chimesOutY.end())-2.0*wrapperHalfXY;
           double yup  = *std::max_element(chimesOutY.begin(),chimesOutY.end())+2.0*wrapperHalfXY;

           G4Box* box = new G4Box("ChimeOut",0.4*outerCaseRadius,wrapperHalfXY,crystalHalfZ);
           G4IntersectionSolid* cutChime1 = new G4IntersectionSolid("Cropped chime", crystalCaseDisk, box,0, G4ThreeVector(0,ylow,0));
           G4IntersectionSolid* cutChime2 = new G4IntersectionSolid("Cropped chime", crystalCaseDisk, box,0, G4ThreeVector(0,yup,0));

           G4LogicalVolume* cutChime1Log = caloBuildLogical(cutChime1, outerRingMaterial, "cutChimeLog",  isDiskVisible, G4Color::Yellow(), isCrystalSolid, forceEdge);   
           G4LogicalVolume* cutChime2Log = caloBuildLogical(cutChime2, outerRingMaterial, "cutChimeLog",  isDiskVisible, G4Color::Yellow(), isCrystalSolid, forceEdge);   

           pv = new G4PVPlacement(0,CLHEP::Hep3Vector(0,0,0),cutChime1Log,"Chimes",crystalCaseDiskLog,true,0,false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
           pv = new G4PVPlacement(0,CLHEP::Hep3Vector(0,0,0),cutChime2Log,"Chimes",crystalCaseDiskLog,true,0,false);
           doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       }
       
       
       
       //------------------------------------------------------------
       // finally put the crystals inside the crystalDiskCase 
       int nTotCrystal(0);
       for (int i=0;i<idisk;++i) nTotCrystal+=cal.disk(idisk).nCrystals();
       
       for (int ic=0; ic <cal.disk(idisk).nCrystals(); ++ic)
       {	      
	     G4int id = nTotCrystal+ic;
	     std::ostringstream name;name<<"caloCrystalPV_" <<id;
             CLHEP::Hep3Vector crystalPosition = cal.disk(idisk).crystal(ic).localPosition();
             crystalPosition.setZ(0.0);

             pv = new G4PVPlacement(0,crystalPosition,wrapperLog,name.str(),crystalCaseDiskLog,true,id,false);
             doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       }

              
       return fullCaseDiskLog;
   }
  
  
  
  
  //--------------------------------------------------------------------------------------------------------------------------------
  // build full backplate - yes this was annoying
  G4LogicalVolume* caloBuildBackPlate(const SimpleConfig& config, MaterialFinder& materialFinder, const DiskCalorimeter& cal, int idisk)
  {            
       const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
       geomOptions->loadEntry( config, "calorimeterRO", "calorimeter.readout" );

       const bool isROVisible          = geomOptions->isVisible("calorimeterRO"); 
       const bool isROSolid            = geomOptions->isSolid("calorimeterRO"); 
       const bool forceEdge             = config.getBool("g4.forceEdge",false);
       const bool doSurfaceCheck        = config.getBool("g4.doSurfaceCheck",false) || config.getBool("calorimeter.doSurfaceCheck",false);
       const int  verbosityLevel        = config.getInt("calorimeter.verbosityLevel",1);

       G4Material* vacuumMaterial       = materialFinder.get("calorimeter.vacuumMaterial");
       G4Material* ROMaterial           = materialFinder.get("calorimeter.readoutMaterial");
       G4Material* FEEMaterial          = materialFinder.get("calorimeter.readoutMaterial");
       G4Material* FEEBoxMaterial       = materialFinder.get("calorimeter.FEEBoxMaterial");
       G4Material* backPlateMaterial    = materialFinder.get("calorimeter.BackPlateMaterial");
       G4Material* pipeMaterial         = materialFinder.get("calorimeter.BackPipeMaterial");
       G4Material* stripMaterial        = materialFinder.get("calorimeter.stripBPMaterial");

       G4double wrapperHalfXY           = cal.caloInfo().crystalXYLength()/2.0 + cal.caloInfo().wrapperThickness();
       G4double ROHalfX                 = cal.caloInfo().roXLength()/2.0;
       G4double ROHalfY                 = cal.caloInfo().roYLength()/2.0;
       G4double ROHalfZ                 = cal.caloInfo().roZLength()/2.0;  
       G4double holeHalfX               = cal.caloInfo().BPHoleXLength()/2.0;
       G4double holeHalfY               = cal.caloInfo().BPHoleYLength()/2.0;
       G4double holeHalfZ               = cal.caloInfo().BPHoleZLength()/2.0;
       G4double FEEHalfX                = cal.caloInfo().FEEXLength()/2.0;
       G4double FEEHalfY                = cal.caloInfo().FEEYLength()/2.0;
       G4double FEEHalfZ                = cal.caloInfo().FEEZLength()/2.0;

       G4double FEEBoxThickness         = cal.caloInfo().FEEBoxThickness();
       G4double FEEBoxHalfX             = holeHalfX + 2*FEEBoxThickness;
       G4double FEEBoxHalfY             = FEEHalfY + 2*FEEBoxThickness;
       G4double FEEBoxHalfZ             = FEEHalfZ + 2*FEEBoxThickness;
       G4double pipeThickness           = cal.caloInfo().pipeThickness();
       G4double coolingPipeRadius       = cal.caloInfo().coolBPPipeRadius();
       G4double stripHalfY              = cal.caloInfo().stripYLength()/2.0;
       G4double stripHalfZ              = cal.caloInfo().stripThickness()/2.0;

       // disk properties
       G4double diskCaseInnerRadius     = cal.geomInfo().diskCaseRadiusIn().at(idisk);
       G4double diskCaseOuterRadius     = cal.geomInfo().diskCaseRadiusOut().at(idisk);
       G4double innerCaseThickness      = cal.caloInfo().caseThicknessIn();
       G4double outerCaseThickness      = cal.caloInfo().caseThicknessOut();
       G4double outerRingEdgeThickness  = cal.caloInfo().outerRingEdgeRLength();     
       G4double innerRingRadius         = diskCaseInnerRadius - innerCaseThickness;
       G4double outerRingRadius         = diskCaseOuterRadius + outerCaseThickness ;   
       G4double outerRingEdgeRadius     = outerRingRadius + outerRingEdgeThickness;   

       std::vector<double> chimesInX   = cal.caloInfo().chimesInsideX();
       std::vector<double> chimesInY   = cal.caloInfo().chimesInsideY();
       std::vector<double> chimesOutX  = cal.caloInfo().chimesOutsideX();
       std::vector<double> chimesOutY  = cal.caloInfo().chimesOutsideY();
       int nChimesInX                  = int(chimesInX.size());
       int nChimesOutX                 = int(chimesOutX.size());


       G4VPhysicalVolume* pv;

       if (cal.caloInfo().nROPerCrystal()==0) return nullptr;


       //-------------------------------------------------------------------------------------------
       // Build hole in back plate, SiPM inside and piece of FEE card (remove to increase speed?)
       //----
       G4Box *holeBack    = new G4Box("caloHoleBack",    holeHalfX,holeHalfY,holeHalfZ);
       G4Box *crystalRO   = new G4Box("caloCrystalRO",   ROHalfX,ROHalfY,ROHalfZ);
       //G4Box *FEEMiniCard = new G4Box("caloFEEMiniCard", FEEHalfX,0.8*holeHalfY,holeHalfZ-ROHalfZ);

       G4LogicalVolume* holeBackLog    = caloBuildLogical(holeBack,   vacuumMaterial, "caloHoleBackLog", isROVisible, G4Color::Black(),isROSolid,forceEdge);   
       G4LogicalVolume* crystalROLog   = caloBuildLogical(crystalRO,  ROMaterial,     "caloCrystalROLog", isROVisible, G4Color::Grey(),isROSolid,forceEdge );        
       //G4LogicalVolume* FEEMiniCardLog = caloBuildLogical(FEEMiniCard,FEEMaterial,    "caloFEEMiniCardLog", isROVisible, G4Color::Green(),isROSolid,forceEdge );        

       // -- place SiPM at bottom of hole    	
       pv = new G4PVPlacement(0,G4ThreeVector(ROHalfX, 0, -holeHalfZ+ROHalfZ),  crystalROLog,"caloROPV_0", holeBackLog, true,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       pv = new G4PVPlacement(0,G4ThreeVector(-ROHalfX, 0, -holeHalfZ+ROHalfZ), crystalROLog,"caloROPV_1", holeBackLog, true,1,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       //pv = new G4PVPlacement(0,G4ThreeVector(ROHalfX, 0, ROHalfZ),  FEEMiniCardLog,"caloFEEMiniCardPV_0", holeBackLog, true,0,false);
       //doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       //pv = new G4PVPlacement(0,G4ThreeVector(-ROHalfX, 0, ROHalfZ), FEEMiniCardLog,"caloFEEMiniCardPV_1", holeBackLog, true,1,false);
       //doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

       // add sensitive detector    
       G4VSensitiveDetector* crSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloReadout());
       if (crSD) crystalROLog->SetSensitiveDetector(crSD);


       //----------------------
       // Build FEE in copper box
       //----
       G4Box *FEEBox    = new G4Box("caloFEEBox",  FEEBoxHalfX,FEEBoxHalfY,FEEBoxHalfZ);
       G4Box *FEEBoxIn  = new G4Box("caloFEEBoxIn",FEEBoxHalfX-FEEBoxThickness,FEEBoxHalfY-FEEBoxThickness,FEEBoxHalfZ-0.5*FEEBoxThickness);
       G4Box *FEECard   = new G4Box("caloFEECard", FEEHalfX,FEEHalfY,FEEHalfZ);

       G4LogicalVolume* FEEBoxLog   = caloBuildLogical(FEEBox,   FEEBoxMaterial, "caloFEEBoxLog",   isROVisible, G4Color::Yellow(),0,forceEdge);   
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

       G4Tubs* backPlate = new G4Tubs("caloBackPlate",innerRingRadius,outerRingRadius,holeHalfZ,0,CLHEP::twopi);       
       G4LogicalVolume* backPlateLog = caloBuildLogical(backPlate, backPlateMaterial, "caloBackPlateLog", isROVisible, G4Color::Blue(),isROSolid,forceEdge);

       for(int ic=0; ic <cal.disk(idisk).nCrystals(); ++ic)
       {	      
	     G4int id = nTotCrystal+ic;
	     std::ostringstream name;name<<"caloHolePV_" <<id;
             CLHEP::Hep3Vector unitPosition = cal.disk(idisk).crystal(ic).localPosition();
             unitPosition.setZ(0.0);

             pv = new G4PVPlacement(0,unitPosition,holeBackLog,name.str(),backPlateLog,false,id,false);
             doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       }

       
             
       // two kind of copper strips: cut by inner hole and full length
       // strips cut by inner hole have return pipes, and adjacent strips need to have the same length
       std::vector<double> stripY, stripX0, stripX1;
       for (int is=0;is<int(outerRingRadius/2.0/wrapperHalfXY);is+=1)
       {
           double ycrystal = is*2.0*wrapperHalfXY;
           
           int idxIn(-1),idxOut(-1);
           for (unsigned i=0; i<chimesInY.size();++i)  if (abs(chimesInY[i]-ycrystal) < 0.1) idxIn=i;       
           for (unsigned i=0; i<chimesOutY.size();++i) if (abs(chimesOutY[i]-ycrystal) < 0.1) idxOut=i;
           
           double x0  = (idxIn > -1) ? chimesInX[idxIn] : 0.0;
           double x1  = (idxOut > -1) ? chimesOutX[idxOut] : 0.0;           
           if (idxIn > -1  && idxIn+1<nChimesInX)   x0 = std::min(x0,chimesInX[idxIn+1]);
           if (idxIn > -1  && idxIn+1==nChimesInX)   x0 = 0.0;
           if (idxOut > -1 && idxOut+1<nChimesOutX) x1 = std::max(x1,chimesOutX[idxOut+1]);

           if (idxIn > -1 && idxIn-1<nChimesInX) std::cout<<is<<" "<<chimesInX[idxIn]<<" "<<chimesInX[idxIn+1]<<" "<<x0<<" "<<chimesInY[idxIn]<<std::endl;
           
           stripY.push_back((2*is+1)*wrapperHalfXY);
           stripX0.push_back(x0);      
           stripX1.push_back(x1);      
       
           //hack to make it look more like the drawing for the current calorimeter - pick one          
           //if (x0>0.1 && is%2==1 ) //if (is<8 && is%2==1 )
           //{
           //   double mi = std::min(stripX0[is],stripX0[is-1]);
           //   double ma = std::max(stripX1[is],stripX1[is-1]);
           //   stripX0[is]  = mi;
           //   stripX0[is-1]= mi;            
           //   stripX1[is]  = ma;
           //   stripX1[is-1]= ma;            
           //}
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
               G4LogicalVolume* bstripLog = caloBuildLogical(bstrip, stripMaterial, "caloStripLog", 1, G4Color::Red(),1,forceEdge);

               pv = new G4PVPlacement(0,G4ThreeVector(0,  stripY[i], zpos), bstripLog,"caloStripPV", backPlateLog, false,0,false);
               doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
  	       pv = new G4PVPlacement(0,G4ThreeVector(0, -stripY[i], zpos), bstripLog,"caloStripPV", backPlateLog, false,0,false);
               doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
           }
       
       }


       //----------------------
       // Build FEE boxes and last cooling pipes
       //----                     
       //NOTE if the design changes and the pipes are farther back, do not forget to change the FEE z position in the loop

       G4Tubs* backPlateFEE = new G4Tubs("caloBackFEEPlate",innerRingRadius,outerRingEdgeRadius,FEEBoxHalfZ,0,CLHEP::twopi);       
       G4LogicalVolume* backPlateFEELog = caloBuildLogical(backPlateFEE, vacuumMaterial, "caloBackPlateFEELog", 0, G4Color::Black(),0,0);

       for(int ic=0; ic <cal.disk(idisk).nCrystals(); ++ic)
       {	      
	     G4int id = nTotCrystal+ic;
	     std::ostringstream name;name<<"caloFEEPV_" <<id;
             CLHEP::Hep3Vector unitPosition = cal.disk(idisk).crystal(ic).localPosition();
             unitPosition.setZ(0); //<--- HERE

             pv = new G4PVPlacement(0,unitPosition,FEEBoxLog,name.str(),backPlateFEELog,false,id,false);
             doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       }

       //just add two big pipes outside and we're done!
       double angMin = std::asin(diskCaseInnerRadius/diskCaseOuterRadius);
       G4Torus* coolBP1 = new G4Torus("coolBP1",coolingPipeRadius-pipeThickness, coolingPipeRadius,outerRingEdgeRadius-3*coolingPipeRadius, 0, 1.5*CLHEP::pi+angMin);
       G4Torus* coolBP2 = new G4Torus("coolBP2",coolingPipeRadius-pipeThickness, coolingPipeRadius,outerRingEdgeRadius-coolingPipeRadius,   0, 1.5*CLHEP::pi+angMin);

       G4LogicalVolume* coolBP1Log = caloBuildLogical(coolBP1, pipeMaterial, "caloCoolBP1Log",isROVisible,G4Color::Blue(),isROSolid,forceEdge);      
       G4LogicalVolume* coolBP2Log = caloBuildLogical(coolBP2, pipeMaterial, "caloCoolBP2Log",isROVisible,G4Color::Red(),isROSolid,forceEdge);

       G4RotationMatrix* rotcoolBP1 = new G4RotationMatrix(CLHEP::HepRotation::IDENTITY);rotcoolBP1->rotateZ(angMin);
       G4RotationMatrix* rotcoolBP2 = new G4RotationMatrix(CLHEP::HepRotation::IDENTITY);rotcoolBP2->rotateZ(0.5*CLHEP::pi);

       pv = new G4PVPlacement(rotcoolBP1,G4ThreeVector(0,0,0),coolBP1Log,"caloCoolBP1PV",backPlateFEELog,false,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
       pv = new G4PVPlacement(rotcoolBP2,G4ThreeVector(0,0,0),coolBP2Log,"caloCoolBP2PV",backPlateFEELog,false,0,false);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);


       //----------------------
       // Make the final full back plate
       //----                     
       G4Tubs* fullBackPlateFEE = new G4Tubs("caloFullBackPlate",innerRingRadius,outerRingEdgeRadius,FEEBoxHalfZ+holeHalfZ,0,CLHEP::twopi);       
       G4LogicalVolume* fullBackPlateFEELog = caloBuildLogical(fullBackPlateFEE, vacuumMaterial, "caloFullBackPlateLog", 0, G4Color::Black(),0,0);

       pv = new G4PVPlacement(0,G4ThreeVector(0,0,-FEEBoxHalfZ),backPlateLog,"caloFullBP1PV",fullBackPlateFEELog,false,0,false);
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
 

       G4double crystalHalfZ             = cal.caloInfo().crystalZLength()/2.0;    
       G4int nBoards                     = cal.caloInfo().nBoard();
       G4double crateDX                  = cal.caloInfo().crateXLength()/2.0;
       G4double crateDY                  = cal.caloInfo().crateYLength()/2.0;
       G4double crateDZ                  = cal.caloInfo().crateZLength()/2.0;
       G4double crateADZ                 = crystalHalfZ;
       G4double crateFShieldDisp         = cal.caloInfo().crateFShieldDeltaZ();
       G4double crateFShieldThick        = cal.caloInfo().crateFShieldThickness();
       G4double crateBottomThick         = cal.caloInfo().crateBShieldThickness();
       G4double crateTopThick            = cal.caloInfo().crateTThickness();
       G4double crateSideThick           = cal.caloInfo().crateSThickness();
       G4double crateFShieldDY           = cal.caloInfo().crateFShieldYLength()/2.0;
       G4double crateFullDZ              = crateDZ+crateFShieldDisp/2.0+crateFShieldThick/2.0;
       G4double crateFullDY              = crateDY+crateBottomThick/2.0;
       G4double crateBoxInDY             = crateDY-crateTopThick/2.0;
       G4double crateBDZ                 = crateFullDZ - crateADZ;

       G4double radiatorDY               = cal.caloInfo().radiatorThickness()/2.0;
       G4double activeStripDY            = cal.caloInfo().activeStripThickness()/2.0;
       G4double passiveStripDY           = cal.caloInfo().passiveStripThickness()/2.0;


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
 
       G4int nCrates            = cal.caloInfo().nCrate();
       G4double phi0Crate       = cal.caloInfo().cratephi0()*CLHEP::degree;
       G4double deltaPhiCrate   = cal.caloInfo().crateDeltaPhi()*CLHEP::degree; 
       G4double crateRadIn      = cal.caloInfo().crateRadiusIn();
       G4double crateRadOut     = crateRadIn+1.005*sqrt(4.0*crate->GetYHalfLength()*crate->GetYHalfLength()+crate->GetXHalfLength()*crate->GetXHalfLength());
       G4double crateHalfLength = crate->GetZHalfLength();   
       G4double cratePosY       = crateRadIn+1.001*crate->GetYHalfLength();    
       G4int nCratesBeforeSpace = cal.caloInfo().nCrateBeforeSpace();

       
       
       
      
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
