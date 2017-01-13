//
// Free function to create the Disk calorimeter.
//
//
// Original author Bertrand Echenard
//
// Notes
//
//  1. a crystal is surrounded by a wrapper, then has radouts and the electronics. Readout and electronics is surrounded by a protective box
//  3. placement (z=0 at the base of the polyhedra, not in the middle of the polyhedra in G4!)
//       - crystal is at z = wrapThickness in wrapping frame
//       - RO is at z = crystalDepth+ROHalfThickness
//  4. Disk holds the wrapped crystals, but the casing is only for the inner / outer radii (no front/back casing yet), so z wrapping = wrapdepth/2.0
//  5. Calibration pipes go in front



#include <array>
#include <iostream>

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/TwoVector.h"

// Mu2e includes.
#include "Mu2eG4/inc/constructDiskCalorimeter.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestTorus.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "Mu2eG4/inc/CaloCrystalSD.hh"
#include "Mu2eG4/inc/CaloReadoutSD.hh"
#include "Mu2eG4/inc/CaloReadoutCardSD.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"

// G4 includes
#include "G4Box.hh"
#include "G4Polyhedra.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"




using namespace std;

namespace mu2e {

  VolumeInfo constructDiskCalorimeter( VolumeInfo const &  mother, SimpleConfig const& config )
  {



    //-- Read parameters from config file
    int const verbosityLevel        = config.getInt("calorimeter.verbosityLevel",0);
    bool const isCalorimeterVisible = config.getBool("calorimeter.calorimeterVisible",false);
    bool const isCalorimeterSolid   = config.getBool("calorimeter.calorimeterSolid",false);
    bool const isDiskCaseVisible    = config.getBool("calorimeter.caseVisible",false);
    bool const isDiskCaseSolid      = config.getBool("calorimeter.caseSolid",false);
    bool const isDiskPipeVisible    = config.getBool("calorimeter.pipeVisible",false);
    bool const isDiskPipeSolid      = config.getBool("calorimeter.pipeSolid",false);
    bool const isCrystalVisible     = config.getBool("calorimeter.crystalVisible",false);
    bool const isCrystalSolid       = config.getBool("calorimeter.crystalSolid",true);
    bool const forceAuxEdgeVisible  = config.getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck       = config.getBool("g4.doSurfaceCheck",false);
    int  const crateVersion         = config.getInt("calorimeter.crateVersion",1);


    //-- A helper class for parsing the config file.
    MaterialFinder materialFinder(config);

    G4Material* diskMaterial    = materialFinder.get("calorimeter.diskMaterial");
    G4Material* fillMaterial    = materialFinder.get("calorimeter.fillMaterial");
    G4Material* crysMaterial    = materialFinder.get("calorimeter.crystalMaterial");
    G4Material* wrapMaterial    = materialFinder.get("calorimeter.crystalWrapper");
    G4Material* readMaterial    = materialFinder.get("calorimeter.readoutMaterial");
    G4Material* pipeMaterial    = materialFinder.get("calorimeter.pipeMaterial");
    G4Material* crateMaterial   = materialFinder.get("calorimeter.crateMaterial");
    G4Material* vacuumMaterial  = materialFinder.get("calorimeter.vacuumMaterial");

    G4Material* ROelectMaterial = materialFinder.get("calorimeter.readoutMaterial");

    G4Material* tempFillerMaterial = materialFinder.get("calorimeter.fillerMaterialName");

    
    G4VPhysicalVolume* pv;

    
    //-- Get all disk /crystal informations here
    //   crystal longitudinal size are defined w.r.t the center of the hexagon (like a radius)
    //   crystal/disk length are defined the usual way (each wrapper adds 2*wrapping_size if all warpped around, wrapping_size otherwise!)
        
    DiskCalorimeter const & cal       = *(GeomHandle<DiskCalorimeter>());

    //calorimeter mother envelope
    G4double mother_inRadius          = cal.caloGeomInfo().envelopeInRadius();
    G4double mother_outRadius         = cal.caloGeomInfo().envelopeOutRadius();
    G4double mother_z0                = cal.caloGeomInfo().envelopeZ0();
    G4double mother_z1                = cal.caloGeomInfo().envelopeZ1();

    // ***************    
    //Temp filler material for particle trapping study
    bool useFiller = config.getBool("calorimeter.useFiller",false);
    double fillerIR = config.getDouble("calorimeter.fillerInnerRadius");
    double fillerOR = config.getDouble("calorimeter.fillerOuterRadius");
    double fillerLen = config.getDouble("calorimeter.fillerFullLength");
    double fillerZC = config.getDouble("calorimeter.fillerZCenter");
    // ***************

    //crystal properties
    G4int    crystalnEdges            = cal.caloGeomInfo().crystalNedges();
    G4double crystalPolysize          = cal.caloGeomInfo().crystalHalfTrans();
    G4double crystalDepth             = 2.0*cal.caloGeomInfo().crystalHalfLength();

    G4double wrapThickness            = cal.caloGeomInfo().wrapperThickness();
    G4double wrapPolysize             = crystalPolysize + wrapThickness;
    G4double wrapDepth                = crystalDepth + wrapThickness; 
    
    
    //readout box properties    
    G4int    nRO                      = cal.caloGeomInfo().nROPerCrystal();
    G4double ROHalfThickness          = cal.caloGeomInfo().roHalfThickness();
    G4double ROHalfTrans              = cal.caloGeomInfo().roHalfTrans();    
    G4double ROElecHalfX              = cal.caloGeomInfo().roElecHalfX();
    G4double ROElecHalfY              = cal.caloGeomInfo().roElecHalfY();
    G4double ROElecHalfZ              = cal.caloGeomInfo().roElecHalfZ();
    G4double ROBoxDepth               = 2*ROHalfThickness +2*ROElecHalfZ;

    
    //crystal + radout properties
    G4double unitDepth                = wrapDepth + ROBoxDepth;
    G4double unitPolysize             = wrapPolysize;
   

    //calorimeter calibration system
    G4int nPipes                      = cal.caloGeomInfo().nPipes();      
    G4double pipeRadius               = cal.caloGeomInfo().pipeRadius();
    G4double pipeThickness            = cal.caloGeomInfo().pipeThickness();
    std::vector<double> pipeTorRadius = cal.caloGeomInfo().pipeTorRadius();


    //calorimeter electronics crates
    G4double crateRadIn               = cal.caloGeomInfo().crateRadiusIn();
    G4double crateRadOut              = cal.caloGeomInfo().crateRadiusOut();
    G4double crateHalfLength          = cal.caloGeomInfo().crateHalfLength();


    //disk properties
    const unsigned int nDisks         = cal.nDisk();
    G4double caseThickness            = cal.caloGeomInfo().caseThickness();
    G4double caseDepth                = unitDepth + 2.0*caseThickness;
    G4double diskDepth                = caseDepth + 2.0*pipeRadius;




    
    
    



    //-----------------------------------------
    // Define and build crystal + wrapper units   
    
	// -- define required solids, need to add an offset to the phi angle to _squares_ to have them rotated properly    
	//
	double offsetAngle = (crystalnEdges==4) ? CLHEP::pi/4.0 : 0; 

	G4double crystalWrapZplanes[2] = {0,wrapDepth};
	G4double crystalWrapRinner[2]  = {0,0};
	G4double crystalWrapRouter[2]  = {wrapPolysize,wrapPolysize};
	G4Polyhedra* crystalWrap       = new G4Polyhedra("CrystalWrap",
                                             offsetAngle,CLHEP::twopi+offsetAngle, 
                                	     crystalnEdges,2, 
					     crystalWrapZplanes,crystalWrapRinner,crystalWrapRouter);

	G4double crystalZplanes[2] = {0,crystalDepth};
	G4double crystalRinner[2]  = {0,0};
	G4double crystalRouter[2]  = {crystalPolysize,crystalPolysize};
	G4Polyhedra* crystal       = new G4Polyhedra("Crystal",
                                	 offsetAngle,CLHEP::twopi+offsetAngle, 
                                	 crystalnEdges,2, 
					 crystalZplanes,crystalRinner,crystalRouter);


	// -- define required logical volumes
	//
	G4LogicalVolume *CrystalLog  = new G4LogicalVolume(crystal, crysMaterial, "CrystalLog");
	G4LogicalVolume *WrapLog = new G4LogicalVolume(crystalWrap, wrapMaterial, "WrapLog");  

	G4VisAttributes* crys_visAtt = new G4VisAttributes(isCrystalVisible, G4Color::Green());
	crys_visAtt->SetForceSolid(isCrystalSolid);
	crys_visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
	CrystalLog->SetVisAttributes(crys_visAtt);    

	WrapLog->SetVisAttributes(G4VisAttributes::Invisible);
	if (isCrystalVisible) WrapLog->SetVisAttributes(G4Color::Magenta());


	// -- place components - refernce point is base of polyhedra -> crystals is wrapperSize away from base
	//
	pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,wrapThickness),CrystalLog,"CrysPV",WrapLog,false,0,false);
	doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);


	// --add sensitive detector    
	//
	G4VSensitiveDetector* ccSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloCrystal());
	if (ccSD) CrystalLog->SetSensitiveDetector(ccSD);



    //------------------------------------------------------------
    // Define and build electronics box with readout + electronics    
       
	// -- define required solids, see comments above
	//
	G4double ROBoxThickness(0);
	
	G4double ROBoxInZplanes[2] = {0,ROBoxDepth-ROBoxThickness};
	G4double ROBoxInRinner[2]  = {0,0};
	G4double ROBoxInRouter[2]  = {crystalPolysize-ROBoxThickness,crystalPolysize-ROBoxThickness};
	G4Polyhedra* ROBoxIn       = new G4Polyhedra("ROBoxIn",
                                       offsetAngle,CLHEP::twopi+offsetAngle, 
                                       crystalnEdges,2, 
                                       ROBoxInZplanes,ROBoxInRinner,ROBoxInRouter);

	/*
	When adding the RO box, make sure to correct the indices for the SD dfor RO and RO Card
	G4double ROBoxOutZplanes[2] = {0,ROBoxDepth};
	G4double ROBoxOutRinner[2]  = {0,0};
	G4double ROBoxOutRouter[2]  = {crystalPolysize,crystalPolysize};
	G4Polyhedra* ROBoxOut       = new G4Polyhedra("ROBoxIn",
                                       offsetAngle,CLHEP::twopi+offsetAngle, 
                                       crystalnEdges,2, 
                                       ROBoxOutZplanes,ROBoxOutRinner,ROBoxOutRouter);
        */
	
	G4Box *crystalRO     = new G4Box("CrystalRO",ROHalfTrans,ROHalfTrans,ROHalfThickness);
	G4Box *electronicsRO = new G4Box("ElectronicsRO",ROElecHalfX,ROElecHalfY,ROElecHalfZ);


	// -- define required logical volumes
	//
	G4LogicalVolume *ROBoxInLog = new G4LogicalVolume(ROBoxIn, vacuumMaterial, "ROBoxInLog");   
	ROBoxInLog->SetVisAttributes(G4VisAttributes::Invisible);
	//if (isCrystalVisible) ROBoxInLog->SetVisAttributes(G4Color::Cyan());

	G4LogicalVolume *ROLog = new G4LogicalVolume(crystalRO, readMaterial, "CrystalROLog" );        
	ROLog->SetVisAttributes(G4VisAttributes::Invisible);
	if (isCrystalVisible) ROLog->SetVisAttributes(G4Color::Grey());

	G4LogicalVolume *ROElectronicsLog = new G4LogicalVolume(electronicsRO, ROelectMaterial, "ElectronicsROLog" );       
	ROElectronicsLog->SetVisAttributes(G4VisAttributes::Invisible);
	if (isCrystalVisible) ROElectronicsLog->SetVisAttributes(G4Color::Green());


	// -- place components - reference point is base of polyhedra and center of RO chip/cards -> RO is on the base+halfsize, 
	//    cards are on top of RO. Card location must be consistent with size    
	//
	double R0disp = 0.5*crystalPolysize;

	std::vector<double> XposRO, YposRO;
	if (nRO==1) {XposRO.push_back(0);
                     YposRO.push_back(0);}
	if (nRO==2) {XposRO.push_back(-R0disp);XposRO.push_back(R0disp); 
                     YposRO.push_back(0);YposRO.push_back(0);}
	if (nRO==4) {XposRO.push_back(-R0disp);XposRO.push_back(-R0disp);XposRO.push_back(R0disp);XposRO.push_back(R0disp);
                     YposRO.push_back(-R0disp);YposRO.push_back(R0disp);YposRO.push_back(-R0disp);YposRO.push_back(R0disp);}

	for (unsigned int iro=0;iro < XposRO.size();++iro)
	{
            ostringstream ROPV; ROPV<<"ROPV_" <<iro;
            pv = new G4PVPlacement(0,G4ThreeVector(XposRO[iro],YposRO[iro],ROHalfThickness), ROLog,ROPV.str(), ROBoxInLog, true,iro,false);
            doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
	}

	//this should be the right position (vertical cards), otherwise invert X-Y and change X-Y size accordingly.	
	pv = new G4PVPlacement(0,G4ThreeVector(R0disp, 0, 2*ROHalfThickness+ROElecHalfZ),  ROElectronicsLog,"ROElectroPV_0", ROBoxInLog, true,0,false);
	doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
	pv = new G4PVPlacement(0,G4ThreeVector(-R0disp, 0, 2*ROHalfThickness+ROElecHalfZ), ROElectronicsLog,"ROElectroPV_1", ROBoxInLog, true,1,false);
	doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);


	// --add sensitive detector    
	//
	G4VSensitiveDetector* crSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloReadout());
	if (crSD) ROLog->SetSensitiveDetector(crSD);
	
	G4VSensitiveDetector* crCardSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloReadoutCard());
	if (crCardSD) ROElectronicsLog->SetSensitiveDetector(crCardSD);
       
       

    //---------------------------------------------------------------
    // Define final crystal + R0box unit here
           
	// -- define required solids
	//
	G4double unitZplanes[2] = {0,unitDepth};
	G4double unitRinner[2]  = {0,0};
	G4double unitRouter[2]  = {unitPolysize,unitPolysize};
	G4Polyhedra* Unit       = new G4Polyhedra("Unit",
                                      offsetAngle,CLHEP::twopi+offsetAngle, 
                                      crystalnEdges,2, 
                                      unitZplanes,unitRinner,unitRouter);


	// -- define required logical volumes
	//
	G4LogicalVolume *UnitLog = new G4LogicalVolume(Unit, diskMaterial, "UnitLog");  //<--CHANGE MATAERIAL 
	UnitLog->SetVisAttributes(G4VisAttributes::Invisible);
	//UnitLog->SetVisAttributes(G4Color::Yellow());

	// -- place components - reference point is base of polyhedra/ wrapper at base, R0Box at WrapDepth
	//
	pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,0.0),WrapLog,"CrysWrapPV",UnitLog,false,0,false);
	doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

	pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,wrapDepth),ROBoxInLog,"ROBoxPV",UnitLog,false,0,false);
	doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

        // Note if we don't add the box around the electronics, we could place directly the RO and the Cards in the Unit, 
	// adding wrapDepth to the Z position of the other guys.








  

    //--------------------------------------
    // Construct calorimeter mother volume
    //
       double mother_zlength  = mother_z1-mother_z0;
       double mother_zCenter  = (mother_z1+mother_z0)/2.0;

       // Make the mother volume for the calorimeter.
       CLHEP::Hep3Vector const& posDS3  = mother.centerInMu2e();
       G4ThreeVector posCaloMother      = G4ThreeVector(posDS3.x(), 0, mother_zCenter);
       G4ThreeVector posCaloMotherInDS  = posCaloMother - posDS3;

       TubsParams caloParams(mother_inRadius,mother_outRadius,mother_zlength/2.0, 0., CLHEP::twopi);
       VolumeInfo calorimeterInfo = nestTubs( "CalorimeterMother",
                                	 caloParams,fillMaterial,0,posCaloMotherInDS,
					 mother,0,
                                	 isCalorimeterVisible,G4Colour::Blue(),isCalorimeterSolid,forceAuxEdgeVisible,
                                	 true,doSurfaceCheck);

       if ( verbosityLevel > 0) 
       {
            double zhl = static_cast<G4Tubs*>(calorimeterInfo.solid)->GetZHalfLength();
	    double CalorimeterOffsetInMu2eZ = calorimeterInfo.centerInMu2e()[CLHEP::Hep3Vector::Z];
	    cout << __func__ << " Calorimeter mother center in Mu2e   : " << calorimeterInfo.centerInMu2e() << endl;
	    cout << __func__ << " Calorimeter mother Z extent in Mu2e : " << CalorimeterOffsetInMu2eZ - zhl << ", " << CalorimeterOffsetInMu2eZ + zhl << endl;
       }



       // ****--------------------------------------------------****
       // Add Temp filler material for particle trapping studies
       if ( useFiller ) {
	 TubsParams fillerParams(fillerIR,fillerOR,fillerLen/2.0);
	 G4ThreeVector posFillerInCaloMother(0.0,0.0,fillerZC-posCaloMother.z());
	 VolumeInfo fillerInfo = nestTubs( "CalorimeterFiller",
                                	 fillerParams,tempFillerMaterial,0,posFillerInCaloMother,
					 calorimeterInfo,0,
                                	 isCalorimeterVisible,G4Colour::Blue(),isCalorimeterSolid,forceAuxEdgeVisible,
                                	 true,doSurfaceCheck);

       }



    //--------------------------------------
    // Construct disks: diskCase contains the crystals. DiskBox contains the calibration pipes and diskCase
    //
       VolumeInfo diskBoxInfo[nDisks];
       VolumeInfo diskCaseInfo[nDisks];
       VolumeInfo diskFEBInfo[nDisks];
       int nTotCrystal(0);

       for (unsigned int idisk=0;idisk<nDisks;++idisk)
       {
	     ostringstream discname0;      discname0<<"DiskCalorimeter_" <<idisk;
	     ostringstream discname1;      discname1<<"DiskCase_" <<idisk;
	     ostringstream cratename;      cratename<<"DiskFEB_" <<idisk;

	     double radiusIn       = cal.disk(idisk).innerRadius()-caseThickness;
	     double radiusOut      = cal.disk(idisk).outerRadius()+caseThickness;
	     double diskpar0[5]    = {radiusIn,radiusOut, diskDepth/2.0, 0., CLHEP::twopi};
	     double diskpar1[5]    = {radiusIn,radiusOut, caseDepth/2.0, 0., CLHEP::twopi};


	     //origin gives the position of the center of the disk, irrespective of the coordinate origin set in the calo description
	     G4ThreeVector posDisk = cal.disk(idisk).origin() - posCaloMother;

	     diskBoxInfo[idisk] =  nestTubs(discname0.str(),
                        		    diskpar0,fillMaterial,&cal.disk(idisk).rotation(),posDisk,
                        		    calorimeterInfo,
                        		    idisk,
                        		    isDiskCaseVisible,G4Colour::Green(),isDiskCaseSolid,forceAuxEdgeVisible,
                        		    true,doSurfaceCheck );

	     diskCaseInfo[idisk] = nestTubs(discname1.str(),
                        		    diskpar1,diskMaterial,0,G4ThreeVector(0.0,0.0,pipeRadius),
                        		    diskBoxInfo[idisk],
                        		    nTotCrystal,
                        		    isDiskCaseVisible,G4Colour::Red(),isDiskCaseSolid,forceAuxEdgeVisible,
                        		    true,doSurfaceCheck );

	     if ( crateVersion > 1 )  // crateVersion 1 is No crates
	     {
		G4ThreeVector posCrate = posDisk + CLHEP::Hep3Vector(0.0,0.0,cal.disk(idisk).crateDeltaZ());
                double cratepar[5] = {crateRadIn, crateRadOut, crateHalfLength, 0., CLHEP::pi};

		diskFEBInfo[idisk] = nestTubs(cratename.str(),
					     cratepar,crateMaterial,&cal.disk(idisk).rotation(),posCrate,
					     calorimeterInfo,
					     idisk,
					     isDiskCaseVisible,G4Colour::Green(),isDiskCaseSolid,forceAuxEdgeVisible,					     
					     true,doSurfaceCheck );
	     }

	     if ( verbosityLevel > 0) 
	     {
		 cout << __func__ << " CalorimeterDisk center in Mu2e    : " << diskBoxInfo[idisk].centerInMu2e() << endl;
		 cout << __func__ << " CalorimeterDisk Z extent in Mu2e  : " << diskBoxInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] - diskDepth/2.0 << ", " << diskBoxInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] + diskDepth/2.0 << endl;
		 cout << __func__ << " CalorimeterCase center in Mu2e    : " << diskCaseInfo[idisk].centerInMu2e() << endl;
		 cout << __func__ << " CalorimeterCase Z extent in Mu2e  : " << diskCaseInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] - caseDepth/2.0 << ", " << diskCaseInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] + caseDepth/2.0 << endl;
	     }


	     // fill this disk with crystal units defined above
	     for(int ic=0; ic <cal.disk(idisk).nCrystals(); ++ic)
	     {	      
        	   CLHEP::Hep3Vector unitPosition = cal.disk(idisk).crystal(ic).localPosition();
        	   double x = unitPosition.x();
        	   double y = unitPosition.y();
        	   double z = -unitDepth/2.0; 	      

		   G4int id = nTotCrystal+ic;
		   ostringstream cryPVName;      
		   cryPVName<<"CrysUnitPV_" <<id;
		   
		   pv = new G4PVPlacement(0,G4ThreeVector(x,y,z),UnitLog,cryPVName.str(),diskCaseInfo[idisk].logical,false,id,false);
        	   doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
	      }

              nTotCrystal += cal.disk(idisk).nCrystals();
	}




     //--------------------------------------
     // Construct Calibration system
     //
	if (nPipes>0 && pipeRadius > 0.001) 
	{

             G4double ZposPipe     = -diskDepth/2.0+pipeRadius;                                            

	     for (unsigned int idisk=0;idisk<nDisks;++idisk)
	     {
		  VolumeInfo caloPipe[nPipes];	
		  for (int ipipe=0;ipipe<nPipes;++ipipe)
		  {
		       ostringstream pipename;  pipename<<"CaloPipe" <<idisk<<"_"<<ipipe;
        	       std::array<double,5> pipeParam { {pipeRadius-pipeThickness, pipeRadius, pipeTorRadius[ipipe], 0., CLHEP::twopi } };

		       caloPipe[ipipe] = nestTorus(pipename.str(),
                                		   pipeParam,pipeMaterial,0,G4ThreeVector(0.0,0.0,ZposPipe),
                                		   diskBoxInfo[idisk],1000*ipipe,
                                		   isDiskPipeVisible,G4Color::Cyan(),isDiskPipeSolid,forceAuxEdgeVisible,
                                		   true,doSurfaceCheck);
		  }				   
             }
	}





     //--------------------------------------
     // All done

     return calorimeterInfo;

  }
  
  

}
