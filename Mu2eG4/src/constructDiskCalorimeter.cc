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
    bool const isDiskBoxVisible     = config.getBool("calorimeter.boxVisible",true);
    bool const isDiskBoxSolid       = config.getBool("calorimeter.boxSolid",true);
    bool const isDiskCaseVisible    = config.getBool("calorimeter.caseVisible",false);
    bool const isDiskCaseSolid      = config.getBool("calorimeter.caseSolid",false);
    bool const isDiskPipeVisible    = config.getBool("calorimeter.pipeVisible",false);
    bool const isDiskPipeSolid      = config.getBool("calorimeter.pipeSolid",false);
    bool const isCrystalVisible     = config.getBool("calorimeter.crystalVisible",false);
    bool const isCrystalSolid       = config.getBool("calorimeter.crystalSolid",true);
    bool const forceAuxEdgeVisible  = config.getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck       = config.getBool("g4.doSurfaceCheck",false);
    bool const placePV              = true;


    //-- A helper class for parsing the config file.
    MaterialFinder materialFinder(config);

    G4Material* diskMaterial   = materialFinder.get("calorimeter.calorimeterDiskMaterial");
    G4Material* fillMaterial   = materialFinder.get("calorimeter.calorimeterFillMaterial");
    G4Material* crysMaterial   = materialFinder.get("calorimeter.crystalMaterial");
    G4Material* wrapMaterial   = materialFinder.get("calorimeter.crystalWrapper");
    G4Material* readMaterial   = materialFinder.get("calorimeter.crystalReadoutMaterial");
    G4Material* pipeMaterial   = materialFinder.get("calorimeter.pipeMaterial");



    
    //-- Get all disk /crystal informations here
    //   crystal longitudinal size are defined w.r.t the center of the hexagon (like a radius)
    //   crystal/disk length are defined the usual way (each wrapper adds 2*wrapping_size!)
        
    DiskCalorimeter const & cal    = *(GeomHandle<DiskCalorimeter>());

    //calorimeter mother enveloppe
    G4double mother_inRadius       = cal.caloGeomInfo().enveloppeInRadius();
    G4double mother_outRadius      = cal.caloGeomInfo().enveloppeOutRadius();
    G4double mother_z0             = cal.caloGeomInfo().enveloppeZ0();
    G4double mother_z1             = cal.caloGeomInfo().enveloppeZ1();
    
    //crystal properties
    G4int    nRO                   = cal.caloGeomInfo().nROPerCrystal();
    G4double ROHalfThickness       = cal.caloGeomInfo().roHalfThickness();
    G4double ROHalfTrans           = cal.caloGeomInfo().roHalfTrans();

    G4double crystalHexsize        = cal.caloGeomInfo().crystalHalfTrans();
    G4double crystalDepth          = 2.0*cal.caloGeomInfo().crystalHalfLength();

    G4double wrapThickness         = cal.caloGeomInfo().wrapperThickness();
    G4double wrapHexsize           = crystalHexsize + wrapThickness;
    G4double wrapDepth             = crystalDepth + 2.0*ROHalfThickness + 2.0*wrapThickness; 
    
    G4double shellThickness        = cal.caloGeomInfo().shellThickness();
    G4double shellHexsize          = wrapHexsize + shellThickness;
    G4double shellDepth            = wrapDepth; //+ 2.0*shellThickness; for a crystal shell surrounding the z faces

    
    // calculate z positions (x,y depend on the crystal)
    G4double crystalZPos           = wrapThickness;
    G4double ROZPos                = wrapThickness+crystalDepth+ROHalfThickness;
    G4double wrapZPos              = (shellDepth-wrapDepth)/2.0;


    //calorimeter calibration system
    G4int nPipes                      = cal.caloGeomInfo().nPipes();      
    G4double pipeRadius               = cal.caloGeomInfo().pipeRadius();
    G4double pipeThickness            = cal.caloGeomInfo().pipeThickness();
    std::vector<double> pipeTorRadius = cal.caloGeomInfo().pipeTorRadius();


    //-- Definition for hexagon dimensions for shell/wrap/crystals   - zplanes, rinner and router define hexagon

    G4double crystalShellZplanes[2] = {0,shellDepth};
    G4double crystalShellRinner[2]  = {0,0};
    G4double crystalShellRouter[2]  = {shellHexsize,shellHexsize};
    G4Polyhedra* crystalShell       = new G4Polyhedra("CrystalShell",
                                          0.,CLHEP::twopi, //phi start phi end
                                          6,2,      //numside ,numzplanes
                                          crystalShellZplanes,crystalShellRinner,crystalShellRouter);

    G4double crystalWrapZplanes[2] = {0,wrapDepth};
    G4double crystalWrapRinner[2]  = {0,0};
    G4double crystalWrapRouter[2]  = {wrapHexsize,wrapHexsize};
    G4Polyhedra* crystalWrap       = new G4Polyhedra("CrystalWrap",
                                	 0.,CLHEP::twopi, //phi start phi end
                                	 6,2,    //numside ,numzplanes
                                	 crystalWrapZplanes,crystalWrapRinner,crystalWrapRouter);

    G4double crystalZplanes[2] = {0,crystalDepth};
    G4double crystalRinner[2]  = {0,0};
    G4double crystalRouter[2]  = {crystalHexsize,crystalHexsize};
    G4Polyhedra* crystal       = new G4Polyhedra("Crystal",
                                     0.,CLHEP::twopi, //phi start phi end
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

    const unsigned int nDisks = cal.nDisk();
    VolumeInfo diskBoxInfo[nDisks];
    VolumeInfo diskCaseInfo[nDisks];
    VolumeInfo diskInnerInfo[nDisks];

    G4VPhysicalVolume* pv;

    //counter of crystals in disks
    G4int crystalIdOffset(0);
    
    for (unsigned int idisk=0;idisk<nDisks;++idisk)
    {

	G4ThreeVector posDisk = cal.origin() + cal.disk(idisk).originLocal() - posCaloMother;

	ostringstream discname0;      discname0<<"DiskCalorimeter_" <<idisk;
	ostringstream discname1;      discname1<<"DiskCase_" <<idisk;
	ostringstream discname2;      discname2<<"DiskInner_" <<idisk;

	double radiusIn   = cal.disk(idisk).innerRadius();
	double radiusOut  = cal.disk(idisk).outerRadius();
	double caseDepth  = shellDepth + 2.0*cal.caloGeomInfo().caseThickness();
	double diskDepth  = caseDepth  + 2.0*pipeRadius;

	double diskpar0[5] = {radiusIn-cal.caloGeomInfo().caseThickness(),radiusOut+cal.caloGeomInfo().caseThickness(), diskDepth/2.0, 0., CLHEP::twopi};
	double diskpar1[5] = {radiusIn-cal.caloGeomInfo().caseThickness(),radiusOut+cal.caloGeomInfo().caseThickness(), caseDepth/2.0, 0., CLHEP::twopi};
	double diskpar2[5] = {radiusIn                    ,radiusOut                    ,shellDepth/2.0, 0., CLHEP::twopi};

	diskBoxInfo[idisk] =  nestTubs(discname0.str(),
                              diskpar0,
                              fillMaterial,
                              &cal.disk(idisk).rotation(),
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
              std::array<double,5> pipeParam { {pipeRadius-pipeThickness, pipeRadius, pipeTorRadius[ipipe], 0., CLHEP::twopi } };

	      caloPipe[ipipe] = nestTorus(pipename.str(),
                                	  pipeParam,
                                	  pipeMaterial,
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


	// -- then fill the inner disk with crystals
	G4int nCrystalInThisDisk = cal.disk(idisk).nCrystals();			
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
              CLHEP::Hep3Vector crystalPosition = cal.disk(idisk).crystal(ic).localPosition();
              double x = crystalPosition.x();
              double y = crystalPosition.y();
              double z = -shellDepth/2.0; 	      

              // place a shell only if it has non-zero thickness, or place the wrapper directly
              if (shellThickness > 0.001)
	      {
	           pv = new G4PVPlacement(0,G4ThreeVector(x,y,z),thisShellLog,"CrysShellPV",
                                          diskInnerInfo[idisk].logical,false,id,false);
                   doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

                   pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,wrapZPos),thisWrapLog,"CrysWrapPV",
                                          thisShellLog,false,id,false);
                   doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
                   
              } else 
	      {
	           pv = new G4PVPlacement(0,G4ThreeVector(x,y,z),thisWrapLog,"CrysWrapPV",
                                          diskInnerInfo[idisk].logical,false,id,false);
                   doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
	      }

	      // -- place crystal inside warp
              pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,crystalZPos),CrystalLog,"CrysPV",
                                     thisWrapLog,false,id,false);
              doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);


              // -- add the readout
	      if (nRO==1) {
                pv = new G4PVPlacement(0,G4ThreeVector(0,0,ROZPos),ROLog,"CrysROPV_0",
                                       thisWrapLog,false,roidBase,false);
                doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
              }
	      if (nRO==2) 
	      { 
		pv = new G4PVPlacement(0,G4ThreeVector(0,-0.5*crystalHexsize,ROZPos),ROLog,"CrysROPV_0",
                                       thisWrapLog,false,roidBase,false);
                doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
        	pv = new G4PVPlacement(0,G4ThreeVector(0, 0.5*crystalHexsize,ROZPos),ROLog,"CrysROPV_1",
                                       thisWrapLog,false,roidBase+1,false);
                doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
              }

	      if (nRO==4) 
	      { 
		pv = new G4PVPlacement(0,G4ThreeVector(-0.5*crystalHexsize,-0.5*crystalHexsize,ROZPos),ROLog,"CrysROPV_0",
                                       thisWrapLog,false,roidBase,false);
                doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

        	pv = new G4PVPlacement(0,G4ThreeVector(-0.5*crystalHexsize, 0.5*crystalHexsize,ROZPos),ROLog,"CrysROPV_1",
                                       thisWrapLog,false,roidBase+1,false);
                doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

        	pv = new G4PVPlacement(0,G4ThreeVector( 0.5*crystalHexsize,-0.5*crystalHexsize,ROZPos),ROLog,"CrysR0PV_2",
                                       thisWrapLog,false,roidBase+2,false);
                doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

        	pv = new G4PVPlacement(0,G4ThreeVector( 0.5*crystalHexsize, 0.5*crystalHexsize,ROZPos),ROLog,"CrysROPV_3",
                                       thisWrapLog,false,roidBase+3,false);
                doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
              }


	   }//end loop over crystals
           
	   crystalIdOffset +=nCrystalInThisDisk;

     }//end loop over disks


     return calorimeterInfo;


  }//end of disk calo construction

} // end namespace mu2e
