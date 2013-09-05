//
// Free function to create the Disk calorimeter.
//
//
// Original author Bertrand Echenard
//
// Notes
//
//  1. a crystal has readouts at the back, both are surrounded by the wrapping, and the wrapping by a shell
//  2. by default, the wrapping surrounds the front/back face of the crystal+ro, the shell does not (shell is a casing)
//  3.  placement (z=0 at the base of the polyhedra, not in the middle of the polyhedra in G4!)
//       - crystal is at z=wrapThickness in wrapping frame
//       - RO is at z=wrapThickness+crystalDepth+ROHalfThickness in wrapping
//       - wrapping is at z=shellThickness in crystal
//  4. Disk holds the crystals but have dead space, so they are enclosed in a larger disk to simulate the casing material

#include <iostream>
#include "CLHEP/Vector/TwoVector.h"

// Mu2e includes.
#include "Mu2eG4/inc/constructHybridCalorimeter.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestTorus.hh"
#include "CalorimeterGeom/inc/HybridCalorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/Barrel.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "Mu2eG4/inc/CaloCrystalSD.hh"
#include "Mu2eG4/inc/CaloReadoutSD.hh"

// G4 includes
#include "G4Box.hh"
#include "G4Polyhedra.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4Trd.hh"



using namespace std;

namespace mu2e {

  VolumeInfo constructHybridCalorimeter( VolumeInfo const &  mother, SimpleConfig const& config )
  {



    //-- Read parameters from config file
    int const verbosityLevel        = config.getInt("calorimeter.verbosityLevel",0);
    bool const isCalorimeterVisible = config.getBool("calorimeter.calorimeterVisible",false);
    bool const isCalorimeterSolid   = config.getBool("calorimeter.calorimeterSolid",false);
    bool const isDiskBoxVisible     = config.getBool("calorimeter.vaneBoxVisible",true);
    bool const isDiskBoxSolid       = config.getBool("calorimeter.vaneBoxSolid",true);
    bool const isDiskCaseVisible    = config.getBool("calorimeter.vaneCaseVisible",false);
    bool const isDiskCaseSolid      = config.getBool("calorimeter.vaneCaseSolid",false);
    bool const isDiskPipeVisible    = config.getBool("calorimeter.vanePipeVisible",false);
    bool const isDiskPipeSolid      = config.getBool("calorimeter.vanePipeSolid",false);
    bool const isCrystalVisible     = config.getBool("calorimeter.crystalVisible",false);
    bool const isCrystalSolid       = config.getBool("calorimeter.crystalSolid",true);
    bool const forceAuxEdgeVisible  = config.getBool("g4.forceAuxEdgeVisible",false);
    G4bool const doSurfaceCheck     = config.getBool("g4.doSurfaceCheck",false);
    bool const placePV              = true;

    //calorimeter mother enveloppe
    double mother_inRadius          = config.getDouble("calorimeter.caloMotherInRadius",0); 
    double mother_outRadius         = config.getDouble("calorimeter.caloMotherOutRadius",765); 
    double mother_z0                = config.getDouble("calorimeter.caloMotherZ0",11740); 
    double mother_z1                = config.getDouble("calorimeter.caloMotherZ1",13910); 

    //calorimeter calibration system
    int const nPipes                = config.getInt("calorimeter.nPipes");      
    double const pipeRadius         = config.getDouble("calorimeter.pipeRadius",5); 
    double const pipeThickness      = config.getDouble("calorimeter.pipeThickness",0.25);     
    std::vector<double> pipeTorRadius;
    config.getVectorDouble("calorimeter.pipeTorRadius",  pipeTorRadius, nPipes);



    //-- A helper class for parsing the config file.
    MaterialFinder materialFinder(config);

    G4Material* diskMaterial   = materialFinder.get("calorimeter.calorimeterDiskMaterial");
    G4Material* fillMaterial   = materialFinder.get("calorimeter.calorimeterFillMaterial");
    G4Material* crysMaterial   = materialFinder.get("calorimeter.crystalMaterial");
    G4Material* wrapMaterial   = materialFinder.get("calorimeter.crystalWrapper");
    G4Material* readMaterial   = materialFinder.get("calorimeter.crystalReadoutMaterial");
    

    
    //-- Get all disk /crystal informations here
    //   crystal longitudinal size are defined w.r.t the center of the hexagon (like a radius)
    //   crystal/disk length are defined the usual way (each wrapper adds 2*wrapping_size!)
        
    HybridCalorimeter const & cal    = *(GeomHandle<HybridCalorimeter>());

    G4int    nRO                   = cal.nROPerCrystal();
    G4double ROHalfThickness       = cal.roHalfThickness();
    G4double ROHalfTrans           = cal.roHalfSize();

    G4double crystalHexsize        = cal.hexCrystalHalfTrans();
    G4double crystalDepth          = 2.0*cal.hexCrystalHalfLength();

    G4double wrapThickness         = cal.wrapperThickness();
    G4double wrapHexsize           = crystalHexsize + wrapThickness;
    G4double wrapDepth             = crystalDepth + 2.0*ROHalfThickness + 2.0*wrapThickness; 
    
    G4double shellThickness        = cal.shellThickness();
    G4double shellHexsize          = wrapHexsize + shellThickness;
    G4double shellDepth            = wrapDepth;// + 2.0*shellThickness;// for a crystal shell surrounding the z faces

    
    // calculate z positions (x,y depend on the crystal)
    G4double crystalZPos           = wrapThickness;
    G4double ROZPos                = wrapThickness+crystalDepth+ROHalfThickness;
    G4double wrapZPos              = (shellDepth-wrapDepth)/2.0;




    //-- Definition for hexagon dimensions for shell/wrap/crystals   - zplanes, rinner and router define hexagon

    G4double crystalShellZplanes[2] = {0,shellDepth};
    G4double crystalShellRinner[2]  = {0,0};
    G4double crystalShellRouter[2]  = {shellHexsize,shellHexsize};
    G4Polyhedra* crystalShell       = new G4Polyhedra("CrystalShell",
						      0,2*pi, //phi start phi end
						      6,2,      //numside ,numzplanes
						      crystalShellZplanes,crystalShellRinner,crystalShellRouter);

    G4double crystalWrapZplanes[2] = {0,wrapDepth};
    G4double crystalWrapRinner[2]  = {0,0};
    G4double crystalWrapRouter[2]  = {wrapHexsize,wrapHexsize};
    G4Polyhedra* crystalWrap       = new G4Polyhedra("CrystalWrap",
						     0,2*pi, //phi start phi end
						     6,2,    //numside ,numzplanes
						     crystalWrapZplanes,crystalWrapRinner,crystalWrapRouter);

    G4double crystalZplanes[2] = {0,crystalDepth};
    G4double crystalRinner[2]  = {0,0};
    G4double crystalRouter[2]  = {crystalHexsize,crystalHexsize};
    G4Polyhedra* crystal       = new G4Polyhedra("Crystal",
						 0,2*pi, //phi start phi end
						 6,2,    //numside//numzplanes
						 crystalZplanes,crystalRinner,crystalRouter);

    //-- Definition of crystal readouts
    G4Box *crystalRO           = new G4Box("CrystalRO",ROHalfTrans,ROHalfTrans,ROHalfThickness);






    //-- Definition of a few logical volumes    
    //
    // Geant4 indexing is such that the crystal /readout can be defined once, 
    // but the shell and wrapper logical volumes must be defined every time a crystal is placed  
    // to get correct index in CaloCrystalSD class

    G4LogicalVolume *CrystalLog  = new G4LogicalVolume(crystal, crysMaterial, "CrystalLog");
    G4LogicalVolume *ROLog       = new G4LogicalVolume(crystalRO, readMaterial, "CrystalROLog" );    

    if(!isCrystalVisible) {
      CrystalLog->SetVisAttributes(G4VisAttributes::Invisible);
      ROLog->SetVisAttributes(G4VisAttributes::Invisible);
    } else {
      G4VisAttributes* crys_visAtt = new G4VisAttributes(isCrystalVisible, G4Color::Green());
      crys_visAtt->SetForceSolid(isCrystalSolid);
      crys_visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
      CrystalLog->SetVisAttributes(crys_visAtt);

      G4VisAttributes* ro_visAtt = new G4VisAttributes(isCrystalVisible, G4Color::Red());
      ro_visAtt->SetForceSolid(isCrystalSolid);
      ro_visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
      ROLog->SetVisAttributes(ro_visAtt);
    }
    
    
    //-- Sensitive detector
    G4VSensitiveDetector* ccSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloCrystal());
    G4VSensitiveDetector* crSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloReadout());
    
    CrystalLog->SetSensitiveDetector(ccSD);
    ROLog->SetSensitiveDetector(crSD);


    //define crystal volume for the Barrel
    G4double tan_barrel_angle  = std::fabs(std::tan(CLHEP::pi / cal.barrel().nCrystalWheel() ) );
    G4double barrel_radiusIn   = cal.barrel().innerRadius();
    G4double barrel_radiusOut  = cal.barrel().outerRadius();
    G4double crystalShell_dx1  = barrel_radiusIn*tan_barrel_angle;// + wrapThickness + shellThickness;
    G4double crystalShell_dy1  = cal.crystalHalfTrans()+ wrapThickness + shellThickness;
    G4double crystalShell_dx2  = barrel_radiusOut*tan_barrel_angle;// + wrapThickness + shellThickness;;
    G4double crystalShell_dy2  = crystalShell_dx1;
    G4double crystalShell_dz   = cal.crystalHalfLength() + ROHalfThickness + wrapThickness + shellThickness;
    G4Trd* barrel_crystalShell = new G4Trd("Barrel_crystalShell",
					   crystalShell_dx1,
					   crystalShell_dx2,
					   crystalShell_dy1,
					   crystalShell_dy2,
					   crystalShell_dz);
    
    G4double crystalWrap_dx1  = barrel_radiusIn*tan_barrel_angle - shellThickness;//+ wrapThickness;
    G4double crystalWrap_dy1  = cal.crystalHalfTrans()+ wrapThickness;
    G4double crystalWrap_dx2  = barrel_radiusOut*tan_barrel_angle - shellThickness;// + wrapThickness;
    G4double crystalWrap_dy2  = crystalWrap_dx1;
    G4double crystalWrap_dz   = cal.crystalHalfLength() + ROHalfThickness + wrapThickness;
    G4Trd* barrel_crystalWrap = new G4Trd("Barrel_crystalWrap",
					  crystalWrap_dx1,
					  crystalWrap_dx2,
					  crystalWrap_dy1,
					  crystalWrap_dy2,
					  crystalWrap_dz);

    G4double crystal_dx1  = barrel_radiusIn*tan_barrel_angle - shellThickness - wrapThickness;//;
    G4double crystal_dy1  = cal.crystalHalfTrans();
    G4double crystal_dx2  = barrel_radiusOut*tan_barrel_angle - shellThickness - wrapThickness;//;
    G4double crystal_dy2  = cal.crystalHalfTrans();
    G4double crystal_dz   = cal.crystalHalfLength();
    G4Trd* barrel_crystal =  new G4Trd("Barrel_crystal",
				       crystal_dx1,
				       crystal_dx2,
				       crystal_dy1,
				       crystal_dy2,
				       crystal_dz);
    
    //-- Definition of crystal readouts
    G4Box *barrel_crystalRO           = new G4Box("Barrel_crystalRO",ROHalfTrans,ROHalfTrans,ROHalfThickness);

    //-- Definition of a few logical volumes    
    //
    // Geant4 indexing is such that the crystal /readout can be defined once, 
    // but the shell and wrapper logical volumes must be defined every time a crystal is placed  
    // to get correct index in CaloCrystalSD class

    G4LogicalVolume *Barrel_crystalLog  = new G4LogicalVolume(barrel_crystal, crysMaterial, "Barrel_crystalLog");
    G4LogicalVolume *Barrel_ROLog       = new G4LogicalVolume(barrel_crystalRO, readMaterial, "Barrel_crystalROLog" );    
    if(!isCrystalVisible) {
      Barrel_crystalLog->SetVisAttributes(G4VisAttributes::Invisible);
      Barrel_ROLog->SetVisAttributes(G4VisAttributes::Invisible);
    } else {
      G4VisAttributes* crys_visAtt = new G4VisAttributes(isCrystalVisible, G4Color::Green());
      crys_visAtt->SetForceSolid(isCrystalSolid);
      crys_visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
      Barrel_crystalLog->SetVisAttributes(crys_visAtt);

      G4VisAttributes* ro_visAtt = new G4VisAttributes(isCrystalVisible, G4Color::Red());
      ro_visAtt->SetForceSolid(isCrystalSolid);
      ro_visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
      Barrel_ROLog->SetVisAttributes(ro_visAtt);
    }
    
    
    //-- Sensitive detector
    // G4VSensitiveDetector* ccSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloCrystal());
    // G4VSensitiveDetector* crSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloReadout());
    
    Barrel_crystalLog->SetSensitiveDetector(ccSD);
    Barrel_ROLog->SetSensitiveDetector(crSD);



    //end description of the barrel components


  

    //-- Construct calorrimeter mother volume
    
    double mother_zlength  = mother_z1-mother_z0;
    double mother_zCenter  = (mother_z1+mother_z0)/2.0;

    //  Make the mother volume for the calorimeter.
    CLHEP::Hep3Vector const& posDS3  = mother.centerInMu2e();
    G4ThreeVector posCaloMother      = G4ThreeVector(posDS3.x(), 0, mother_zCenter);
    G4ThreeVector posCaloMotherInDS  = posCaloMother - posDS3;


    TubsParams caloParams(mother_inRadius,mother_outRadius,mother_zlength/2.0, 0., CLHEP::twopi);
    VolumeInfo calorimeterInfo = nestTubs( "CalorimeterMother",
					   caloParams,
					   fillMaterial,
					   0,
					   posCaloMotherInDS,
					   mother,
					   0,
					   isCalorimeterVisible,
					   G4Colour::Blue(),
					   isCalorimeterSolid,
					   forceAuxEdgeVisible,
					   placePV,
					   doSurfaceCheck
					   );

    if ( verbosityLevel > 0) 
      {
	double zhl         = static_cast<G4Tubs*>(calorimeterInfo.solid)->GetZHalfLength();
	CLHEP::Hep3Vector const & CalorimeterOffsetInMu2e = calorimeterInfo.centerInMu2e();
	double CalorimeterOffsetInMu2eZ = CalorimeterOffsetInMu2e[CLHEP::Hep3Vector::Z];
	cout << __func__ << " Calorimeter mother center in Mu2e   : " << CalorimeterOffsetInMu2e << endl;
	cout << __func__ << " Calorimeter mother Z extent in Mu2e : " << CalorimeterOffsetInMu2eZ - zhl << ", " << CalorimeterOffsetInMu2eZ + zhl << endl;
      }



    //-- Construct disks: diskInner contains the crystals, and is inside diskCase. DiskBox contains the calibration pipes and diskCase

    const unsigned int nDisks = cal.nSections();
    VolumeInfo diskBoxInfo[nDisks];
    VolumeInfo diskCaseInfo[nDisks];
    VolumeInfo diskInnerInfo[nDisks];


    //counter of crystals in disks
    G4int crystalIdOffset(0);
    int idisk(0);
    

    G4ThreeVector posDisk = cal.origin() + cal.disk().originLocal() - posCaloMother;

    ostringstream discname0;      discname0<<"HybridCalorimeter_";
    ostringstream discname1;      discname1<<"HybridCase_";
    ostringstream discname2;      discname2<<"HybridInner_";

    double radiusIn   = cal.disk().innerRadius();
    double radiusOut  = cal.disk().outerRadius();
    double caseDepth  = shellDepth + 2.0*cal.caseThickness();
    double diskDepth  = caseDepth  + 2.0*pipeRadius;

    double diskpar0[5] = {radiusIn-cal.caseThickness(),radiusOut+cal.caseThickness(), diskDepth/2.0, 0, 2*pi};
    double diskpar1[5] = {radiusIn-cal.caseThickness(),radiusOut+cal.caseThickness(), caseDepth/2.0, 0, 2*pi};
    double diskpar2[5] = {radiusIn                    ,radiusOut                    ,shellDepth/2.0, 0, 2*pi};

    diskBoxInfo[idisk] =  nestTubs(discname0.str(),
				   diskpar0,
				   fillMaterial,
				   &cal.disk().rotation(),
				   posDisk,
				   calorimeterInfo,
				   idisk,
				   isDiskCaseVisible,
				   G4Colour::Green(),
				   isDiskCaseSolid,
				   forceAuxEdgeVisible,
				   placePV,
				   doSurfaceCheck );

    diskCaseInfo[idisk] = nestTubs(discname1.str(),
				   diskpar1,
				   diskMaterial,
				   0,
				   G4ThreeVector(0.0,0.0,pipeRadius),
				   diskBoxInfo[idisk],
				   10*idisk,
				   isDiskCaseVisible,
				   G4Colour::Red(),
				   isDiskCaseSolid,
				   forceAuxEdgeVisible,
				   placePV,
				   doSurfaceCheck );

    diskInnerInfo[idisk] = nestTubs(discname2.str(),
				    diskpar2,
				    diskMaterial,
				    0,
				    G4ThreeVector(0.0,0.0,0.0),
				    diskCaseInfo[idisk],
				    100*idisk,
				    isDiskBoxVisible,
				    G4Colour::Green(),
				    isDiskBoxSolid,
				    forceAuxEdgeVisible,
				    placePV,
				    doSurfaceCheck );
			      

    if ( verbosityLevel > 0) 
      {
	cout << __func__ << " CalorimeterDisk center in Mu2e    : " << diskBoxInfo[idisk].centerInMu2e() << endl;
	cout << __func__ << " CalorimeterDisk Z extent in Mu2e  : " << diskBoxInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] - diskDepth/2.0 << ", " << diskBoxInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] + diskDepth/2.0 << endl;
	cout << __func__ << " CalorimeterCase center in Mu2e    : " << diskCaseInfo[idisk].centerInMu2e() << endl;
	cout << __func__ << " CalorimeterCase Z extent in Mu2e  : " << diskCaseInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] - caseDepth/2.0 << ", " << diskCaseInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] + caseDepth/2.0 << endl;
	cout << __func__ << " CalorimeterInner center in Mu2e    : " << diskInnerInfo[idisk].centerInMu2e() << endl;
	cout << __func__ << " CalorimeterInner Z extent in Mu2e  : " << diskInnerInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] - shellDepth/2.0 << ", " << diskInnerInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] + shellDepth/2.0 << endl;
      }



    //add pipes from calibration system
    if (nPipes>0 && pipeRadius>0.01) 
      {
	   VolumeInfo caloPipe[nPipes];	
	   for (int ipipe=0;ipipe<nPipes;++ipipe){

	      ostringstream pipename;  pipename<<"CaloPipe" <<idisk<<"_"<<ipipe;
              std::array<double,5> pipeParam { {pipeRadius-pipeThickness, pipeRadius, pipeTorRadius[ipipe], 0, 2*pi } };

	      caloPipe[ipipe] = nestTorus(pipename.str(),
                                	  pipeParam,
                                	  diskMaterial,
                                	  0,
                                	  G4ThreeVector(0.0,0.0,-diskDepth/2.0+pipeRadius),
                                	  diskBoxInfo[idisk],
                                	  1000*ipipe,
                                	  isDiskPipeVisible,
                                	  G4Color::Cyan(),
                                	  isDiskPipeSolid,
                                	  forceAuxEdgeVisible,
                                	  placePV,
                                	  doSurfaceCheck
                                	  );
	   }				   
        }

    if ( verbosityLevel > 0) 
      {
	printf("\nPie volume created\n");
	printf("start filling crystals\n");
	
      }
    // -- then fill the inner disk with crystals
    G4int nCrystalInThisDisk = cal.disk().nCrystals();			
    for(int ic=0; ic <nCrystalInThisDisk; ++ic)
      {

	G4int id       = crystalIdOffset+ic;
	G4int roidBase = cal.ROBaseByCrystal(id);

	// Have to define a shell / wrapper logical volume for each crystal 
	// to get correct index in CrystalCaloSD
	G4LogicalVolume *thisShellLog(0);
	if (shellThickness > 0.001)
	  {
	    thisShellLog = new G4LogicalVolume(crystalShell, fillMaterial, "ShellLog");
	    thisShellLog->SetVisAttributes(G4VisAttributes::Invisible);
	  }

	G4LogicalVolume *thisWrapLog = new G4LogicalVolume(crystalWrap, wrapMaterial, "WrapLog");
	thisWrapLog->SetVisAttributes(G4VisAttributes::Invisible);

	//position of shell in the disk
	//contrary to rectangles, z position of hexagon is their base, not their center in Geant 4!!	      
	CLHEP::Hep3Vector crystalPosition = cal.disk().crystal(ic).position();
	double x = crystalPosition.x();
	double y = crystalPosition.y();
	double z = -shellDepth/2.0; 	      

	// place a shell only if it has non-zero thickness, or place the wrapper directly
	if (shellThickness > 0.001)
	  {
	    new G4PVPlacement(0,G4ThreeVector(x,y,z),thisShellLog,"CrysShellPV",diskInnerInfo[idisk].logical,0,id,doSurfaceCheck);   
	    new G4PVPlacement(0,G4ThreeVector(0.0,0.0,wrapZPos),thisWrapLog,"CrysWrapPV",thisShellLog,0,id,doSurfaceCheck);
	  } else 
	  {
	    new G4PVPlacement(0,G4ThreeVector(x,y,z),thisWrapLog,"CrysWrapPV",diskInnerInfo[idisk].logical,0,id,doSurfaceCheck);   	      
	  }

	// -- place crystal inside warp
	new G4PVPlacement(0,G4ThreeVector(0.0,0.0,crystalZPos),CrystalLog,"CrysPV",thisWrapLog,0,id,doSurfaceCheck);


	// -- add the readout
	if (nRO==1) 
	  new G4PVPlacement(0,G4ThreeVector(0,0,ROZPos),ROLog,"CrysROPV_0",thisWrapLog,0,roidBase,doSurfaceCheck);

	if (nRO==2) 
	  { 
	    new G4PVPlacement(0,G4ThreeVector(0,-0.5*crystalHexsize,ROZPos),ROLog,"CrysROPV_0",thisWrapLog,0,roidBase,doSurfaceCheck);
	    new G4PVPlacement(0,G4ThreeVector(0, 0.5*crystalHexsize,ROZPos),ROLog,"CrysROPV_1",thisWrapLog,0,roidBase+1,doSurfaceCheck);
	  }

	if (nRO==4) 
	  { 
	    new G4PVPlacement(0,G4ThreeVector(-0.5*crystalHexsize,-0.5*crystalHexsize,ROZPos),ROLog,"CrysROPV_0",thisWrapLog,0,roidBase,doSurfaceCheck);
	    new G4PVPlacement(0,G4ThreeVector(-0.5*crystalHexsize, 0.5*crystalHexsize,ROZPos),ROLog,"CrysROPV_1",thisWrapLog,0,roidBase+1,doSurfaceCheck);
	    new G4PVPlacement(0,G4ThreeVector( 0.5*crystalHexsize,-0.5*crystalHexsize,ROZPos),ROLog,"CrysR0PV_2",thisWrapLog,0,roidBase+2,doSurfaceCheck);
	    new G4PVPlacement(0,G4ThreeVector( 0.5*crystalHexsize, 0.5*crystalHexsize,ROZPos),ROLog,"CrysROPV_3",thisWrapLog,0,roidBase+3,doSurfaceCheck);
	  }


      }//end loop over crystals
           
    if ( verbosityLevel > 0) 
      {
	printf("\nend disk calorimeter crystal filling\n");
	printf("start barrel creation\n");
      }
    crystalIdOffset +=nCrystalInThisDisk;

        
    ++idisk;
   
    G4ThreeVector posBarrel = cal.origin() + cal.barrel().originLocal() - posCaloMother;

    ostringstream barrelname0;      barrelname0<<"HybridCalorimeter_" <<idisk;
    ostringstream barrelname1;      barrelname1<<"HybridCase_" <<idisk;
    ostringstream barrelname2;      barrelname2<<"HybridInner_" <<idisk;
  
    G4int nCrystalInThisBarrel = cal.barrel().nCrystals();	
    double Angle_step = CLHEP::twopi / cal.barrel().nCrystalWheel();
    double tmpAngle(0);
    
    radiusIn   = cal.barrel().innerRadius();
    radiusOut  = cal.barrel().outerRadius()*(1.0 + std::tan(Angle_step*0.5)*std::sin(Angle_step*0.5)) ;
    
    double barrelWheels = cal.barrel().nWheels();
    
    double barrelCaseDepth  = barrelWheels*2.0*crystalShell_dy1 + 2.0*cal.caseThickness();//shellDepth + 2.0*cal.caseThickness();
    double barrelShellDepth = barrelWheels*2.0*crystalShell_dy1;
    //double barrelDepth  = caseDepth  + 2.0*pipeRadius;
    
    if ( verbosityLevel > 0) 
      {
	printf("\nBarrel info:\nradiusIn = %.1f,\nradiusOut = %.1f, \nbarelCaseDepth = %.1f, \nbarrelShellDepth = %.1f\n",radiusIn, radiusOut, barrelCaseDepth, barrelShellDepth);
      }
    double dist = cal.caseThickness();// + wrapThickness + shellThickness + ROHalfThickness*2.0 ;
    double barrelpar0[5] = {radiusIn,//cal.caseThickness(),
			    radiusOut,// + dist ,// cal.caseThickness(),
			    barrelCaseDepth/2.0,
			    0, 2*pi};
    if ( verbosityLevel > 0) 
      {
	printf("barrelpar0 = {%.1f, %.1f, %.1f, %.1f, %.1f}\n",
	       barrelpar0[0], barrelpar0[1], barrelpar0[2], barrelpar0[3],
	       barrelpar0[4]);
      }
    double barrelpar1[5] = {radiusIn, //cal.caseThickness(),
			    radiusOut,// + dist,//cal.caseThickness(),
			    barrelCaseDepth/2.0,
			    0, 2*pi};
    if ( verbosityLevel > 0) 
      {
	printf("barrelpar1 = {%.1f, %.1f, %.1f, %.1f, %.1f}\n",
	       barrelpar1[0], barrelpar1[1], barrelpar1[2], barrelpar1[3],
	       barrelpar1[4]);
      }
    double barrelpar2[5] = {radiusIn,
			    radiusOut - dist,
			    barrelShellDepth/2.0,
			    0, 2*pi};
    
    if ( verbosityLevel > 0) 
      {
	printf("barrelpar1 = {%.1f, %.1f, %.1f, %.1f, %.1f}\n",
	       barrelpar2[0], barrelpar2[1], barrelpar2[2], barrelpar2[3],
	       barrelpar2[4]);
	std::cout<<"creating barrelBoxInfo"<<std::endl;
      }
    diskBoxInfo[idisk] =  nestTubs(barrelname0.str(),
				   barrelpar0,
				   fillMaterial,
				   &cal.barrel().rotation(),
				   posBarrel,
				   calorimeterInfo,
				   idisk,
				   isDiskCaseVisible,
				   G4Colour::Green(),
				   isDiskCaseSolid,
				   forceAuxEdgeVisible,
				   placePV,
				   doSurfaceCheck );
    if ( verbosityLevel > 0) 
      {
	printf("\ncreating barrelCaseInfo");
      }
    diskCaseInfo[idisk] = nestTubs(barrelname1.str(),
				   barrelpar1,
				   diskMaterial,
				   0,
				   G4ThreeVector(0.0,0.0,pipeRadius),
				   diskBoxInfo[idisk],
				   10*idisk,
				   isDiskCaseVisible,
				   G4Colour::Red(),
				   isDiskCaseSolid,
				   forceAuxEdgeVisible,
				   placePV,
				   doSurfaceCheck );
    if ( verbosityLevel > 0) 
      {
	printf("\ncreating barrelInnerInfo");
      }
    diskInnerInfo[idisk] = nestTubs(barrelname2.str(),
				    barrelpar2,
				    diskMaterial,
				    0,
				    G4ThreeVector(0.0,0.0,0.0),
				    diskCaseInfo[idisk],
				    100*idisk,
				    isDiskBoxVisible,
				    G4Colour::Green(),
				    isDiskBoxSolid,
				    forceAuxEdgeVisible,
				    placePV,
				    doSurfaceCheck );
			      

    if ( verbosityLevel > 0) 
      {
	cout << __func__ << " BarrelCalorimeter center in Mu2e    : " << diskBoxInfo[idisk].centerInMu2e() << endl;
	cout << __func__ << " BarrelCalorimeter Z extent in Mu2e  : " << diskBoxInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] - crystal_dy1 << ", " << diskBoxInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] + crystal_dy1 << endl;
	cout << __func__ << " BarrelCalorimeterCase center in Mu2e    : " << diskCaseInfo[idisk].centerInMu2e() << endl;
	cout << __func__ << " BarrelCalorimeterCase Z extent in Mu2e  : " << diskCaseInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] - barrelCaseDepth/2.0 << ", " << diskCaseInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] + barrelCaseDepth/2.0 << endl;
	cout << __func__ << " BarrelCalorimeterInner center in Mu2e    : " << diskInnerInfo[idisk].centerInMu2e() << endl;
	cout << __func__ << " BarrelCalorimeterInner Z extent in Mu2e  : " << diskInnerInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] - barrelShellDepth/2.0 << ", " << diskInnerInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] + barrelShellDepth/2.0 << endl;
      }



    //FIX ME
    //add pipes from calibration system
    //      if (nPipes>0 && pipeRadius>0.01) 
    // {
    // 	  VolumeInfo caloPipe[nPipes];	
    // 	  for (int ipipe=0;ipipe<nPipes;++ipipe){

    // 	    ostringstream pipename;  pipename<<"CaloPipe" <<idisk<<"_"<<ipipe;
    // 	    double pipeParam[5] = {pipeRadius-pipeThickness, pipeRadius, pipeTorRadius[ipipe], 0, 2*pi };

    // 	    caloPipe[ipipe] = nestTorus(pipename.str(),
    // 					pipeParam,
    // 					diskMaterial,
    // 					0,
    // 					G4ThreeVector(0.0,0.0,-diskDepth/2.0+pipeRadius),
    // 					diskBoxInfo[idisk],
    // 					1000*ipipe,
    // 					isDiskPipeVisible,
    // 					G4Color::Cyan(),
    // 					isDiskPipeSolid,
    // 					forceAuxEdgeVisible,
    // 					placePV,
    // 					doSurfaceCheck
    // 					);
    // 	  }				   
    //         }


  
    
   
    for(int ic=0; ic < nCrystalInThisBarrel; ++ic)
      {

	G4int id       = crystalIdOffset+ic;
	G4int roidBase = cal.ROBaseByCrystal(id);

	// Have to define a shell / wrapper logical volume for each crystal 
	// to get correct index in CrystalCaloSD
	G4LogicalVolume *thisShellLog(0);
	if (shellThickness > 0.001)
	  {
	    thisShellLog = new G4LogicalVolume(barrel_crystalShell, fillMaterial, "ShellLog");
	    thisShellLog->SetVisAttributes(G4VisAttributes::Invisible);
	  }

	G4LogicalVolume *thisWrapLog = new G4LogicalVolume(barrel_crystalWrap, wrapMaterial, "WrapLog");
	thisWrapLog->SetVisAttributes(G4VisAttributes::Invisible);

	//position of shell in the disk
	//contrary to rectangles, z position of hexagon is their base, not their center in Geant 4!!	      
	int wheelIndex = int(ic/cal.barrel().nCrystalWheel());
	CLHEP::Hep3Vector crystalPosition = cal.barrel().crystal(ic).position();
	double x = crystalPosition.x();
	double y = crystalPosition.y();
	double z = -( barrelShellDepth/2.0 - crystalShell_dy1*(1.0 + 2.0*double(wheelIndex)) );


	//printf("\nFILLING CRYSTALS");
	//printf("\ncrystal position: (%.1f, %.1f, %.1f)", x, y, z);
	//DEFINE THE ROTATION
	int cry_module = ic % int(cal.barrel().nCrystalWheel());
	tmpAngle = double(cry_module)*Angle_step;
	G4RotationMatrix *rotCry = new G4RotationMatrix();//RForTrapezoids);

	rotCry->rotateZ(tmpAngle);
	rotCry->rotateX(90*deg);

	// place a shell only if it has non-zero thickness, or place the wrapper directly
	if (shellThickness > 0.001)
	  {
	    new G4PVPlacement(rotCry
			      , G4ThreeVector(x,y,z)
			      ,thisShellLog
			      ,"CrysShellPV"
			      ,diskInnerInfo[idisk].logical
			      ,0
			      ,id
			      ,doSurfaceCheck);   
	    new G4PVPlacement(0
			      ,G4ThreeVector(0.0,0.0,wrapZPos)/*G4ThreeVector(wrapZPos*std::cos(tmpAngle), wrapZPos*std::sin(tmpAngle), 0.0)*/
			      ,thisWrapLog
			      ,"CrysWrapPV"
			      ,thisShellLog
			      ,0
			      ,id
			      ,doSurfaceCheck);
	  } else 
	  {
	    new G4PVPlacement(rotCry
			      ,G4ThreeVector(x,y,z)
			      ,thisWrapLog
			      ,"CrysWrapPV"
			      ,diskInnerInfo[idisk].logical
			      ,0
			      ,id
			      ,doSurfaceCheck);   	      
	  }

	new G4PVPlacement(0
 			  ,G4ThreeVector(0.0,0.0,crystalZPos)
			  ,Barrel_crystalLog
 			  ,"CrysPV"
 			  ,thisWrapLog
 			  ,0
 			  ,id
			  ,doSurfaceCheck);


	// -- add the readout
	ROZPos = crystalWrap_dz;// + ROHalfThickness;//wrapThickness+crystal_dz+ROHalfThickness;
	//ROZPos *= -1.;
	if (nRO==1) 
	  new G4PVPlacement(0
			    ,G4ThreeVector(0,0,ROZPos)
			    ,Barrel_ROLog
			    ,"CrysROPV_0"
			    ,thisWrapLog
			    ,0
			    ,roidBase
			    ,doSurfaceCheck);

	if (nRO==2) 
	  { 
	    new G4PVPlacement(0, G4ThreeVector(-0.5*crystal_dx2, 0.0, ROZPos),Barrel_ROLog,"CrysROPV_0",thisWrapLog,0,roidBase,doSurfaceCheck);
	    new G4PVPlacement(0, G4ThreeVector( 0.5*crystal_dx2, 0., ROZPos),Barrel_ROLog,"CrysROPV_1",thisWrapLog,0,roidBase+1,doSurfaceCheck);
	  }

	if (nRO==4) 
	  { 
	    new G4PVPlacement(0, G4ThreeVector(-0.5*crystal_dx2,-0.5*crystal_dy2,ROZPos),Barrel_ROLog,"CrysROPV_0",thisWrapLog,0,roidBase,doSurfaceCheck);
	    new G4PVPlacement(0, G4ThreeVector(-0.5*crystal_dx2, 0.5*crystal_dy2,ROZPos),Barrel_ROLog,"CrysROPV_1",thisWrapLog,0,roidBase+1,doSurfaceCheck);
	    new G4PVPlacement(0, G4ThreeVector( 0.5*crystal_dx2,-0.5*crystal_dy2,ROZPos),Barrel_ROLog,"CrysR0PV_2",thisWrapLog,0,roidBase+2,doSurfaceCheck);
	    new G4PVPlacement(0, G4ThreeVector( 0.5*crystal_dx2, 0.5*crystal_dy2,ROZPos),Barrel_ROLog,"CrysROPV_3",thisWrapLog,0,roidBase+3,doSurfaceCheck);
	  }


      }//end loop over crystals
           
    //crystalIdOffset +=nCrystalInThisDisk;

    return calorimeterInfo;


  }//end of disk calo construction

} // end namespace mu2e
