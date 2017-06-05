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

#include "Mu2eG4/inc/constructDiskCalorimeter.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
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

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/TwoVector.h"

#include <array>
#include <iostream>
#include <sstream>



namespace mu2e {

  VolumeInfo constructDiskCalorimeter( const VolumeInfo&  mother, const SimpleConfig& config )
  {



    //-- Read parameters from config file
    int const verbosityLevel        = config.getInt("calorimeter.verbosityLevel",1);
    bool const isCalorimeterVisible = config.getBool("calorimeter.calorimeterVisible",false);
    bool const isCalorimeterSolid   = config.getBool("calorimeter.calorimeterSolid",false);
    bool const isDiskCaseVisible    = config.getBool("calorimeter.caseVisible",false);
    bool const isDiskCaseSolid      = config.getBool("calorimeter.caseSolid",false);
    bool const isDiskPipeVisible    = config.getBool("calorimeter.pipeVisible",false);
    bool const isDiskPipeSolid      = config.getBool("calorimeter.pipeSolid",false);
    bool const isCrystalVisible     = config.getBool("calorimeter.crystalVisible",false);
    bool const isCrystalSolid       = config.getBool("calorimeter.crystalSolid",true);
    bool const forceAuxEdgeVisible  = config.getBool("g4.forceAuxEdgeVisible",false);
    bool const isShieldSolid        = config.getBool("calorimeter.shieldSolid",false);
    bool const doSurfaceCheck       = config.getBool("g4.doSurfaceCheck",false) || config.getBool("calorimeter.doSurfaceCheck",false);
    int  const crateVersion         = config.getInt("calorimeter.crateVersion",1);


    //-- A helper class for parsing the config file.
    MaterialFinder materialFinder(config);

    G4Material* diskMaterial          = materialFinder.get("calorimeter.diskMaterial");
    G4Material* fillMaterial          = materialFinder.get("calorimeter.fillMaterial");
    G4Material* crysMaterial          = materialFinder.get("calorimeter.crystalMaterial");
    G4Material* wrapMaterial          = materialFinder.get("calorimeter.crystalWrapper");
    G4Material* readMaterial          = materialFinder.get("calorimeter.readoutMaterial");
    G4Material* pipeMaterial          = materialFinder.get("calorimeter.pipeMaterial");
    G4Material* crateMaterial         = materialFinder.get("calorimeter.crateMaterial");
    G4Material* crateBkgMaterial      = materialFinder.get("calorimeter.crateBackgroundMaterial");
    G4Material* crateBottomAMaterial  = materialFinder.get("calorimeter.crateMaterial");
    G4Material* shieldMaterial        = materialFinder.get("calorimeter.shieldMaterial");
    // G4Material* copperPlaneMaterial   = materialFinder.get("calorimeter.copperPlaneMaterial");
    // G4Material* substrateMaterial     = materialFinder.get("calorimeter.substrateMaterial");
    G4Material* radiatorMaterial      = materialFinder.get("calorimeter.radiatorMaterial");
    G4Material* activeStripMaterial   = materialFinder.get("calorimeter.activeStripMaterial");
    G4Material* passiveStripMaterial  = materialFinder.get("calorimeter.passiveStripMaterial");
    G4Material* innerRingMaterial     = materialFinder.get("calorimeter.innerRingMaterial");
    G4Material* outerRingMaterial     = materialFinder.get("calorimeter.outerRingMaterial");
    G4Material* outerRingEdgeMaterial = materialFinder.get("calorimeter.outerRingEdgeMaterial");
    // G4Material* stepsMaterial         = materialFinder.get("calorimeter.stepsMaterial");
    G4Material* vacuumMaterial        = materialFinder.get("calorimeter.vacuumMaterial");

    G4Material* ROelectMaterial = materialFinder.get("calorimeter.readoutMaterial");
    G4Material* tempFillerMaterial = materialFinder.get("calorimeter.fillerMaterialName");
    

    G4VPhysicalVolume* pv;

    
    //-- Get all disk /crystal informations here
    //   crystal longitudinal size are defined w.r.t the center of the hexagon (like a radius)
    //   crystal/disk length are defined the usual way (each wrapper adds 2*wrapping_size if all warpped around, wrapping_size otherwise!)
        
    DiskCalorimeter const & cal       = *(GeomHandle<DiskCalorimeter>());

    //calorimeter mother envelope
    G4double mother_inRadius          = cal.caloInfo().envelopeInRadius();
    G4double mother_outRadius         = cal.caloInfo().envelopeOutRadius();
    G4double mother_z0                = cal.caloInfo().envelopeZ0();
    G4double mother_z1                = cal.caloInfo().envelopeZ1();

    // ***************    
    //Temp filler material for particle trapping study
    bool useFiller = config.getBool("calorimeter.useFiller",false);
    double fillerIR = config.getDouble("calorimeter.fillerInnerRadius");
    double fillerOR = config.getDouble("calorimeter.fillerOuterRadius");
    double fillerLen = config.getDouble("calorimeter.fillerFullLength");
    double fillerZC = config.getDouble("calorimeter.fillerZCenter");
    // ***************

    //crystal properties
    G4int    crystalnEdges            = 4;
    G4double crystalPolysize          = cal.caloInfo().crystalHalfTrans();
    G4double crystalDepth             = 2.0*cal.caloInfo().crystalHalfLength();

    G4double wrapThickness            = cal.caloInfo().wrapperThickness();
    G4double wrapPolysize             = crystalPolysize + wrapThickness;
    G4double wrapDepth                = crystalDepth + wrapThickness; 
    
    
    //readout box properties    
    G4int    nRO                      = cal.caloInfo().nROPerCrystal();
    G4double ROHalfThickness          = cal.caloInfo().roHalfThickness();
    G4double ROHalfTrans              = cal.caloInfo().roHalfTrans();    
    G4double ROElecHalfX              = cal.caloInfo().roElecHalfX();
    G4double ROElecHalfY              = cal.caloInfo().roElecHalfY();
    G4double ROElecHalfZ              = cal.caloInfo().roElecHalfZ();
    G4double ROBoxDepth               = 2*ROHalfThickness +2*ROElecHalfZ;

    
    //crystal + radout properties
    G4double unitDepth                = wrapDepth + ROBoxDepth;
    G4double unitPolysize             = wrapPolysize;

    //calorimeter calibration system
    G4int nPipes                      = cal.caloInfo().nPipes();      
    G4double pipeRadius               = cal.caloInfo().pipeRadius();
    G4double pipeThickness            = cal.caloInfo().pipeThickness();
    std::vector<double> pipeTorRadius = cal.caloInfo().pipeTorRadius();
    
    // disk properties
    const unsigned int nDisks         = cal.nDisk();
    G4double stepsInnerRadius         = cal.caloInfo().stepsRadiusIn();
    G4double stepsOuterRadius         = cal.caloInfo().stepsRadiusOut();
    G4double caseThickness            = cal.caloInfo().caseThickness();
    G4double caseInnerThickness       = cal.caloInfo().caseThicknessIn();
    G4double caseOuterThickness       = cal.caloInfo().caseThicknessOut();
    G4double outerRingEdgeDepth       = cal.caloInfo().outerRingEdgeDepth();
    G4double outerRingEdgeThickness   = cal.caloInfo().outerRingEdgeThickness(); 
    
    G4double caseDepth                = unitDepth + 2.0*caseThickness;
    G4double diskDepth                = caseDepth + 2.0*pipeRadius;
    G4double ringDepth                = crystalDepth - 2.0*wrapThickness;
    G4double innerRingRadius          = stepsInnerRadius - caseInnerThickness;
    G4double outerRingRadius          = stepsOuterRadius + caseOuterThickness;   
    G4double outerRingEdgeRadius      = outerRingRadius + outerRingEdgeThickness;
    G4double delta                    = 0.01;


    //z position of the inner and the outer ring
    G4double ZposRing     = pipeRadius - ROBoxDepth/2.0;
    G4double ZposEdge     = pipeRadius + (outerRingEdgeDepth - unitDepth)/2.0;
    
    // crate properties
    G4int nCrates      = cal.nCrate();
    G4int nBoards      = cal.nBoard();

    G4double crateDX   = cal.caloInfo().crateHalfX();
    G4double crateDY   = cal.caloInfo().crateHalfY();
    G4double crateDZ   = cal.caloInfo().crateHalfZ();
    G4double crateADZ  = wrapDepth/2.;
    G4double crateToDiskDeltaZ = cal.disk(0).geomInfo().crateDeltaZ();
    G4double crateBDZ  = (crateDZ + crateToDiskDeltaZ - crateADZ)/2.;

    G4double crateDYTop    = cal.caloInfo().crateTopHalfY();
    G4double crateDXSide   = cal.caloInfo().crateSideHalfX();
    G4double crateDYSide   = cal.caloInfo().crateSideHalfY()-delta;
    G4double crateDYBottom = cal.caloInfo().crateBottomHalfY();

    G4double crateRadIn      = cal.caloInfo().crateRadiusIn();
    G4double crateRadOut     = cal.caloInfo().crateRadiusOut();
    G4double crateHalfLength = cal.caloInfo().crateHalfLength();   

    // shielding properties   
    G4double crateFShieldDistanceDZ  = 23.5;
    G4double shieldDispZ       = cal.caloInfo().shieldDispZ();//60.;//26.5; 
    G4double shieldDZFront     = cal.caloInfo().shieldDZFront(); 
    
    G4double shieldDYFront    = cal.caloInfo().shieldHalfY();
    //    G4double shieldPosYBottom = -crateDY;
    G4double shieldPosYFront  = -crateDY+delta  + shieldDYFront;
    
        
    //    G4double shieldDZBottom   = crateDZ - shieldDZFront;   
    G4double shieldPosZFront  = -crateDZ - shieldDispZ + crateFShieldDistanceDZ + 2.*shieldDZFront;// - shieldDZBottom;
    
    
    // crate coordinates
    G4double crateTopPosY    = crateDY-crateDYTop+delta;
    G4double crateTopPosZ    = crateFShieldDistanceDZ;      
    G4double crateSidePosX   = crateDX-crateDXSide;
    G4double crateSidePosY   = crateDY-2*crateDYTop-crateDYSide+delta;
    G4double crateSidePosZ   = crateFShieldDistanceDZ;    
    G4double crateBottomPosY = crateDYBottom-crateDY+delta;
    //    G4double crateBottomPosZ = crateFShieldDistanceDZ;
    G4double crateBottomAPosZ = -crateToDiskDeltaZ + crateFShieldDistanceDZ;
    G4double crateBottomBPosZ = crateBottomAPosZ + crateADZ + crateBDZ;

    // boards properties
    G4double boardDY        = cal.caloInfo().boardHalfY();
    G4double radiatorDY     = cal.caloInfo().radiatorHalfY();
    G4double activeStripDY  = cal.caloInfo().activeStripHalfY();
    G4double passiveStripDY = cal.caloInfo().passiveStripHalfY();
    
    G4double boardDX    = crateDX-2.*crateDXSide-2*delta;    
    G4double boardDZ    = crateDZ;
    G4double boardDispY = 2*(crateDY-crateDYTop-crateDYBottom)/nBoards;
    G4double boardPosY  = crateBottomPosY+crateDYBottom+boardDispY/2.;
    G4double boardPosZ  = crateFShieldDistanceDZ;    
    // G4double radiatorPosY = boardDY-radiatorDY;   
    // G4double activeStripPosY  = radiatorPosY-radiatorDY-activeStripDY;
    // G4double passiveStripPosY = activeStripPosY - activeStripDY - passiveStripDY;
    G4double radiatorPosY = -boardDY + radiatorDY;   
    G4double activeStripPosY  = radiatorPosY+radiatorDY+activeStripDY;
    G4double passiveStripPosY = activeStripPosY + activeStripDY + passiveStripDY;

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

	G4VisAttributes* crys_visAtt = new G4VisAttributes(isCrystalVisible, G4Color::Black());
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
	// if (isCrystalVisible) ROBoxInLog->SetVisAttributes(G4Color::Cyan());

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
            std::ostringstream ROPV; ROPV<<"ROPV_" <<iro;
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
	    std::cout << __func__ << " Calorimeter mother center in Mu2e   : " << calorimeterInfo.centerInMu2e() << std::endl;
	    std::cout << __func__ << " Calorimeter mother Z extent in Mu2e : " << CalorimeterOffsetInMu2eZ - zhl << ", " << CalorimeterOffsetInMu2eZ + zhl << std::endl;
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
       VolumeInfo calorimeterFEBInfo[nDisks];

       VolumeInfo diskBoxInfo[nDisks];
       VolumeInfo diskCaseInfo[nDisks];
       VolumeInfo diskInnerRingInfo[nDisks];
       VolumeInfo diskOuterRingInfo[nDisks];
       // VolumeInfo diskStepsInfo[nDisks];
       VolumeInfo diskOuterRingFstEdgeInfo[nDisks];
       VolumeInfo diskOuterRingSndEdgeInfo[nDisks];
       int nTotCrystal(0);

       for (unsigned int idisk=0;idisk<nDisks;++idisk)
       {
	     std::ostringstream discname0;      discname0<<"DiskCalorimeter_" <<idisk;
	     std::ostringstream discname1;      discname1<<"DiskCase_" <<idisk;
	     std::ostringstream discname2;      discname2<<"DiskOuterRing_" <<idisk;
	     // std::ostringstream discname3;      discname3<<"DiskSteps_" <<idisk;
	     std::ostringstream discname4;      discname4<<"DiskInnerRing_" <<idisk;
	     std::ostringstream discname5;      discname5<<"DiskOuterRingFstEdge_" <<idisk;
	     std::ostringstream discname6;      discname6<<"DiskOuterRingSndEdge_" <<idisk;

	     std::ostringstream cratename;      cratename<<"CalorimeterFEB_" <<idisk;

	     double radiusIn       = cal.disk(idisk).innerRadius()-caseThickness;
	     double radiusOut      = cal.disk(idisk).outerRadius()+caseThickness;
	     
	     double diskpar0[5]    = {innerRingRadius,outerRingEdgeRadius, diskDepth/2.0, 0., CLHEP::twopi};
	     double diskpar1[5]    = {radiusIn,radiusOut, caseDepth/2.0, 0., CLHEP::twopi};
	     double diskpar2[5]    = {innerRingRadius,stepsInnerRadius, ringDepth/2.0, 0., CLHEP::twopi};
	     // double diskpar3[5]    = {stepsInnerRadius+delta,stepsOuterRadius, ringDepth/2.0, 0., CLHEP::twopi};
	     double diskpar4[5]    = {stepsOuterRadius,outerRingRadius, ringDepth/2.0, 0., CLHEP::twopi};
	     double diskpar5[5]    = {outerRingRadius+delta,outerRingEdgeRadius, outerRingEdgeDepth/2.0, 0., CLHEP::twopi};
	     double diskpar6[5]    = {outerRingRadius+delta,outerRingEdgeRadius, outerRingEdgeDepth/2.0, 0., CLHEP::twopi};

	     //origin gives the position of the center of the disk, irrespective of the coordinate origin set in the calo description
	     G4ThreeVector posDisk = cal.disk(idisk).geomInfo().origin() - posCaloMother;

	     diskBoxInfo[idisk] =  nestTubs(discname0.str(),
                        		    diskpar0,fillMaterial,&cal.disk(idisk).geomInfo().rotation(),posDisk,
                        		    calorimeterInfo,
                        		    idisk,
                        		    isDiskCaseVisible,G4Colour::Cyan(),isDiskCaseSolid,forceAuxEdgeVisible,
                        		    true,doSurfaceCheck );

	     diskCaseInfo[idisk] = nestTubs(discname1.str(),
                        		    diskpar1,diskMaterial,0,G4ThreeVector(0.0,0.0,pipeRadius),
                        		    diskBoxInfo[idisk],
                        		    nTotCrystal,
                        		    isDiskCaseVisible,G4Colour::Red(),isDiskCaseSolid,forceAuxEdgeVisible,
                        		    true,doSurfaceCheck );

	     diskOuterRingInfo[idisk] = nestTubs(discname2.str(),
						 diskpar2,outerRingMaterial,0,G4ThreeVector(0.0,0.0,ZposRing),
						 diskBoxInfo[idisk] ,
						 idisk,
						 isDiskCaseVisible,G4Colour::Blue(),isDiskCaseSolid,forceAuxEdgeVisible,
						 true,doSurfaceCheck );

	     // diskStepsInfo[idisk] = nestTubs(discname3.str(),
	     // 				     diskpar3,stepsMaterial,0,G4ThreeVector(0.0,0.0,ZposRing),
	     // 				     diskBoxInfo[idisk],
	     // 				     idisk,
	     // 				     isDiskCaseVisible,G4Colour::Black(),isDiskCaseSolid,forceAuxEdgeVisible,
	     // 				     true,doSurfaceCheck );
	     
	     diskInnerRingInfo[idisk] = nestTubs(discname4.str(),
						 diskpar4,innerRingMaterial,0,G4ThreeVector(0.0,0.0,ZposRing),
						 diskBoxInfo[idisk],
						 idisk,
						 isDiskCaseVisible,G4Colour::Blue(),isDiskCaseSolid,forceAuxEdgeVisible,
						 true,doSurfaceCheck );
	     
	     diskOuterRingFstEdgeInfo[idisk] = nestTubs(discname5.str(),
	     						diskpar5,outerRingEdgeMaterial,0,G4ThreeVector(0.0,0.0,ZposEdge),
	     						diskBoxInfo[idisk],
	     						idisk,
	     						isDiskCaseVisible,G4Colour::Cyan(),isDiskCaseSolid,forceAuxEdgeVisible,
	     						true,doSurfaceCheck);

	     diskOuterRingSndEdgeInfo[idisk] = nestTubs(discname6.str(),
	     						diskpar6,outerRingEdgeMaterial,0,G4ThreeVector(0.0,0.0,-ZposEdge-(outerRingEdgeDepth+ROBoxDepth)/2.0),
	     						 diskBoxInfo[idisk],
	     						 idisk,
	     						 isDiskCaseVisible,G4Colour::Cyan(),isDiskCaseSolid,forceAuxEdgeVisible,
	     						 true,doSurfaceCheck);

	     if ( crateVersion > 1 )  // crateVersion 1 is No crates
	     {
	       // define the crate box
	       G4Box *crateBox     = new G4Box("CrateBox",crateDX+delta,crateDY+delta,crateDZ+shieldDispZ/2.+shieldDZFront+delta);
	       G4LogicalVolume *crateBoxLog = new G4LogicalVolume(crateBox, vacuumMaterial, "crateBoxLog");   
	       crateBoxLog->SetVisAttributes(G4Color::Black());

	       // define the crate walls
	       G4Box *crateTop     = new G4Box("CrateTop",crateDX,crateDYTop,crateDZ);
	       G4Box *crateSide    = new G4Box("CrateSide",crateDXSide,crateDYSide,crateDZ);
	       G4Box *crateBottomA = new G4Box("CrateBottomA",crateDX-delta,crateDYBottom,crateADZ);
	       G4Box *crateBottomB = new G4Box("CrateBottomB",crateDX-delta,crateDYBottom,crateBDZ);
	       
	       G4LogicalVolume *crateTopLog = new G4LogicalVolume(crateTop, crateMaterial, "crateTopLog");   
	       crateTopLog->SetVisAttributes(G4Color::Black());
	       pv = new G4PVPlacement(0,G4ThreeVector(0.0,crateTopPosY,crateTopPosZ),crateTopLog,"crateTopPV",crateBoxLog,false,0,false);
	       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
	       
	       G4LogicalVolume *crateSideLog = new G4LogicalVolume(crateSide, crateMaterial, "crateSideLog");   
	       crateSideLog->SetVisAttributes(G4Color::Black());
	       pv = new G4PVPlacement(0,G4ThreeVector(crateSidePosX,crateSidePosY,crateSidePosZ),crateSideLog,"crateSidePV_0",crateBoxLog,false,0,false);
	       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
	       pv = new G4PVPlacement(0,G4ThreeVector(-crateSidePosX,crateSidePosY,crateSidePosZ),crateSideLog,"crateSidePV_1",crateBoxLog,false,1,false);
	       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);

	       if (isShieldSolid)
		 {
		   crateBottomAMaterial = shieldMaterial;
		 }
	       else 
		 {
		   crateBottomAMaterial = crateMaterial;
		 }
	       G4LogicalVolume *crateBottomALog = new G4LogicalVolume(crateBottomA, crateBottomAMaterial, "crateBottomALog");   
	       crateBottomALog->SetVisAttributes(G4Color::Black());
	       pv = new G4PVPlacement(0,G4ThreeVector(0.0,crateBottomPosY,crateBottomAPosZ),crateBottomALog,"crateBottomAPV",crateBoxLog,false,0,false);
	       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);

	       G4LogicalVolume *crateBottomBLog = new G4LogicalVolume(crateBottomB, crateMaterial, "crateBottomBLog");   
	       crateBottomBLog->SetVisAttributes(G4Color::Black());
	       pv = new G4PVPlacement(0,G4ThreeVector(0.0,crateBottomPosY,crateBottomBPosZ),crateBottomBLog,"crateBottomBPV",crateBoxLog,false,0,false);
	       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);

	       // define boards
	       G4Box *boardCrate = new G4Box("BoardCrate",boardDX,boardDY,boardDZ);
	       G4LogicalVolume *boardCrateLog = new G4LogicalVolume(boardCrate, vacuumMaterial, "boardCrateLog");   
	       boardCrateLog->SetVisAttributes(G4Color::Black());
	       
	       // define board slices components
	       G4Box *radiatorBoard     = new G4Box("RadiatorBoard",boardDX,radiatorDY,boardDZ);
	       G4Box *activeStripBoard  = new G4Box("ActiveStripBoard",boardDX,activeStripDY,boardDZ);
	       G4Box *passiveStripBoard = new G4Box("PassiveStripBoard",boardDX,passiveStripDY,boardDZ);

	       G4LogicalVolume *radiatorBoardLog = new G4LogicalVolume(radiatorBoard, radiatorMaterial, "radiatorBoardLog");   
	       radiatorBoardLog->SetVisAttributes(G4Color::Black());
	       pv = new G4PVPlacement(0,G4ThreeVector(0.0,radiatorPosY,0.0),radiatorBoardLog,"radiatorBoardPV",boardCrateLog,false,0,false);
	       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);

	       G4LogicalVolume *activeStripBoardLog = new G4LogicalVolume(activeStripBoard, activeStripMaterial, "activeStripBoardLog");   
	       activeStripBoardLog->SetVisAttributes(G4Color::Black());
	       pv = new G4PVPlacement(0,G4ThreeVector(0.0,activeStripPosY,0.0),activeStripBoardLog,"activeStripBoardPV",boardCrateLog,false,0,false);
	       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);

	       G4LogicalVolume *passiveStripBoardLog = new G4LogicalVolume(passiveStripBoard, passiveStripMaterial, "passiveStripBoardLog");   
	       passiveStripBoardLog->SetVisAttributes(G4Color::Black());
	       pv = new G4PVPlacement(0,G4ThreeVector(0.0,passiveStripPosY,0.0),passiveStripBoardLog,"passiveStripBoardPV",boardCrateLog,false,0,false);
	       doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);

	       // put boards in a single crate
	       boardPosY  = crateBottomPosY+crateDYBottom+boardDispY/2.;

	       for (G4int ibrd=0;ibrd < nBoards;++ibrd)
		 {
		   std::ostringstream boardPV; boardPV<<"boardPV_" <<ibrd;

		   G4VSensitiveDetector* crCrate = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloCrate());
		   if (crCrate) activeStripBoardLog->SetSensitiveDetector(crCrate);

		   pv = new G4PVPlacement(0,G4ThreeVector(0.0,boardPosY,boardPosZ), boardCrateLog, boardPV.str(), crateBoxLog, true, ibrd,false);
		   doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
		   boardPosY += boardDispY;
		 }

	       if (isShieldSolid)
		 {
		   // define shielding
		   G4Box           *shieldCrateFront    = new G4Box("shieldCrateFront",crateDX,shieldDYFront,shieldDZFront);
		   G4LogicalVolume *shieldCrateFrontLog = new G4LogicalVolume(shieldCrateFront, shieldMaterial, "ShieldCrateTopLog");   
		   shieldCrateFrontLog->SetVisAttributes(G4Color::Black());
		   
		   pv = new G4PVPlacement(0,G4ThreeVector(0.0,shieldPosYFront,shieldPosZFront),shieldCrateFrontLog,"shieldCrateBottomPV",crateBoxLog,false,0,false);
		   doSurfaceCheck && checkForOverlaps(pv,config,verbosityLevel>0);
		 }		  
	       
	       // put crates onto the calorimeter disks
	       G4double cratePosY     = crateRadIn+crateDY+2.*crateDYBottom + 4*delta;//20*delta;
	       G4double phi0Crate     = (30./360.)*2*CLHEP::pi;
	       G4double deltaPhiCrate = (16./360.)*2*CLHEP::pi;
	       G4double phiCrate(0);

		G4int numberOfCratesBeforeFSpace = 3;

		G4ThreeVector cratePosZ = posDisk + CLHEP::Hep3Vector(0., 0., crateHalfLength - crateADZ);//CLHEP::Hep3Vector(0.0,0.0, crateToDiskDeltaZ );
                double cratepar[5] = {crateRadIn, crateRadOut, crateHalfLength+10*delta,
				      -phi0Crate, CLHEP::pi+2*phi0Crate};

		calorimeterFEBInfo[idisk] = nestTubs(cratename.str(),
						     cratepar,
						     crateBkgMaterial,
						     &cal.disk(idisk).geomInfo().rotation(),posCaloMotherInDS+cratePosZ,
						     mother,
						     idisk,
						     isDiskCaseVisible,
						     G4Colour::Green(),
						     isDiskCaseSolid,forceAuxEdgeVisible,true,doSurfaceCheck);		

		// Dave Brown (Lou) invading code here to put in cable runs
		if ( config.getBool("ds.hasCableRunCal",false)) {

		  double crRin = config.getDouble("ds.CableRunCal.Rin")*CLHEP::mm;
		  double crRout = config.getDouble("ds.CableRunCal.Rout")*CLHEP::mm;

		  std::ostringstream crTubName;
		  crTubName << "CableRunCalTub" << idisk;
		  G4Tubs* ccrTub = new G4Tubs( crTubName.str(),crRin, crRout,
					       crateHalfLength - 5.0,
					       config.getDouble("ds.CableRunCal.phi0")*CLHEP::degree,
					       config.getDouble("ds.CableRunCal.dPhi")*CLHEP::degree);


		  CLHEP::Hep3Vector calCableRunLoc(0.0,0.0,0.0);
		  std::ostringstream ccrLogName;
		  ccrLogName << "CalCableRunLogInCalFeb" << idisk;
		
		  G4LogicalVolume* ccrTubLog = 
		    new G4LogicalVolume(ccrTub,findMaterialOrThrow(config.getString("ds.CableRunCal.material")),ccrLogName.str());

		  std::ostringstream ccrName;
		  ccrName << "CalCableRunInCalFeb" << idisk;
		  G4RotationMatrix ccrRot = G4RotationMatrix();
		  G4Transform3D ccrCoord = G4Transform3D(ccrRot,calCableRunLoc);

		  pv = new G4PVPlacement(ccrCoord,
					 ccrTubLog,ccrName.str(),
					 calorimeterFEBInfo[idisk].logical, 
					 false, 0, false);
		   doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
			     
		} // End of conditionally adding CalCableRun

		if ( config.getBool("ds.hasCableRunTrk",false)) {

		  double crRin = config.getDouble("ds.CableRunTrk.Rin")*CLHEP::mm;
		  double crRout = config.getDouble("ds.CableRunTrk.Rout")*CLHEP::mm;
		  double phi01 = config.getDouble("ds.CableRunTrk.phi0")*CLHEP::degree;
		  double dPhi = config.getDouble("ds.CableRunTrk.dPhi")*CLHEP::degree;
		  double phi02 = 180.0*CLHEP::degree - phi01 - dPhi;
		  if ( phi01 + dPhi > 180.0*CLHEP::degree + phi0Crate ) {
		    dPhi = 179.5*CLHEP::degree + phi0Crate - phi01;
		    phi02 = 180.0*CLHEP::degree - phi01 - dPhi;
		  }
		    
		  std::ostringstream cr1TubName;
		  cr1TubName << "CableRun1TrkTub" << idisk;
		  G4Tubs* ccr1Tub = new G4Tubs( cr1TubName.str(),crRin, crRout,
					       crateHalfLength - 5.0,
						phi01, dPhi);

		  std::ostringstream cr2TubName;
		  cr2TubName << "CableRun2TrkTub" << idisk;
		  G4Tubs* ccr2Tub = new G4Tubs( cr2TubName.str(),crRin, crRout,
					       crateHalfLength - 5.0,
						phi02, dPhi);



		  CLHEP::Hep3Vector trkCableRunLoc(0.0,0.0,0.0);
		  std::ostringstream ccr1LogName;
		  ccr1LogName << "TrkCableRun1LogInCalFeb" << idisk;
		  std::ostringstream ccr2LogName;
		  ccr2LogName << "TrkCableRun2LogInCalFeb" << idisk;
		
		  G4LogicalVolume* ccr1TubLog = 
		    new G4LogicalVolume(ccr1Tub,findMaterialOrThrow(config.getString("ds.CableRunTrk.material")),ccr1LogName.str());
		  G4LogicalVolume* ccr2TubLog = 
		    new G4LogicalVolume(ccr2Tub,findMaterialOrThrow(config.getString("ds.CableRunTrk.material")),ccr2LogName.str());

		  std::ostringstream ccr1Name;
		  ccr1Name << "TrkCableRun1InCalFeb" << idisk;
		  std::ostringstream ccr2Name;
		  ccr2Name << "TrkCableRun2InCalFeb" << idisk;

		  G4RotationMatrix ccrRot = G4RotationMatrix();
		  G4Transform3D ccrCoord = G4Transform3D(ccrRot,trkCableRunLoc);

		  pv = new G4PVPlacement(ccrCoord,
					 ccr1TubLog,ccr1Name.str(),
					 calorimeterFEBInfo[idisk].logical, 
					 false, 0, false);
		   doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

		  pv = new G4PVPlacement(ccrCoord,
					 ccr2TubLog,ccr2Name.str(),
					 calorimeterFEBInfo[idisk].logical, 
					 false, 0, false);
		   doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
			     
		} // End of conditionally adding TrkCableRun


		G4RotationMatrix rotCrate = G4RotationMatrix();
		
		for (G4int icrt=0;icrt < nCrates;++icrt)
		 {		   
		   
		   if (icrt<nCrates/2) 
		     {
		       if (icrt == 0) phiCrate = CLHEP::pi/2. - deltaPhiCrate;	      
		       else 
			 {
			   phiCrate -= deltaPhiCrate;
			   if (icrt == numberOfCratesBeforeFSpace) phiCrate -= deltaPhiCrate;
			 }
		     }
		   else
		     {
		       if (icrt == nCrates/2) phiCrate = CLHEP::pi/2. + deltaPhiCrate;		       
		       else
			 {
			   phiCrate += deltaPhiCrate;
			   if (icrt == nCrates/2+numberOfCratesBeforeFSpace) phiCrate += deltaPhiCrate;
			 }
		     }
		   		   
		   if (icrt==nCrates/2) rotCrate.rotateZ(CLHEP::pi/2.+deltaPhiCrate/3.);
 
		   if (icrt<nCrates/2) 
		     {
		       rotCrate.rotateZ(-deltaPhiCrate);
		       if (icrt == numberOfCratesBeforeFSpace) rotCrate.rotateZ(-deltaPhiCrate);
		     }
		   else 
		     {
		       rotCrate.rotateZ(deltaPhiCrate);
		       if (icrt == nCrates/2+numberOfCratesBeforeFSpace) rotCrate.rotateZ(deltaPhiCrate);
		     }

		   G4ThreeVector posCrate   = G4ThreeVector(cratePosY*std::cos(phiCrate),cratePosY*std::sin(phiCrate),0.0);
		   G4Transform3D crateCoord = G4Transform3D(rotCrate,posCrate);
		   
		   std::ostringstream cratePV; cratePV<<"cratePV_" <<icrt;
		   pv = new G4PVPlacement(crateCoord, crateBoxLog, cratePV.str(), calorimeterFEBInfo[idisk].logical, true, icrt,false);
		   doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
		 }
	     }

	     if ( verbosityLevel > 0) 
	     {
		 std::cout << __func__ << " CalorimeterDisk center in Mu2e    : " << diskBoxInfo[idisk].centerInMu2e() << std::endl;
		 std::cout << __func__ << " CalorimeterDisk Z extent in Mu2e  : " << diskBoxInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] - diskDepth/2.0 << ", " << diskBoxInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] + diskDepth/2.0 << std::endl;
		 std::cout << __func__ << " CalorimeterCase center in Mu2e    : " << diskCaseInfo[idisk].centerInMu2e() << std::endl;
		 std::cout << __func__ << " CalorimeterCase Z extent in Mu2e  : " << diskCaseInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] - caseDepth/2.0 << ", " << diskCaseInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] + caseDepth/2.0 << std::endl;
	     }


	     // fill this disk with crystal units defined above
	     for(int ic=0; ic <cal.disk(idisk).nCrystals(); ++ic)
	     {	      
        	   CLHEP::Hep3Vector unitPosition = cal.disk(idisk).crystal(ic).localPosition();
        	   double x = unitPosition.x();
        	   double y = unitPosition.y();
        	   double z = -unitDepth/2.0; 	      

		   G4int id = nTotCrystal+ic;
		   std::ostringstream cryPVName;      
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
		       std::ostringstream pipename;  pipename<<"CaloPipe" <<idisk<<"_"<<ipipe;
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
