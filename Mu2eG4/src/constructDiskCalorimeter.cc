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
    bool const isDiskCaseVisible    = config.getBool("calorimeter.caseVisible",false);
    bool const isDiskCaseSolid      = config.getBool("calorimeter.caseSolid",false);
    bool const isDiskPipeVisible    = config.getBool("calorimeter.pipeVisible",false);
    bool const isDiskPipeSolid      = config.getBool("calorimeter.pipeSolid",false);
    bool const isCrystalVisible     = config.getBool("calorimeter.crystalVisible",false);
    bool const isCrystalSolid       = config.getBool("calorimeter.crystalSolid",true);
    bool const forceAuxEdgeVisible  = config.getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck       = config.getBool("g4.doSurfaceCheck",false);


    //-- A helper class for parsing the config file.
    MaterialFinder materialFinder(config);

    G4Material* diskMaterial   = materialFinder.get("calorimeter.calorimeterDiskMaterial");
    G4Material* fillMaterial   = materialFinder.get("calorimeter.calorimeterFillMaterial");
    G4Material* crysMaterial   = materialFinder.get("calorimeter.crystalMaterial");
    G4Material* wrapMaterial   = materialFinder.get("calorimeter.crystalWrapper");
    G4Material* readMaterial   = materialFinder.get("calorimeter.crystalReadoutMaterial");
    G4Material* pipeMaterial   = materialFinder.get("calorimeter.pipeMaterial");

    G4VPhysicalVolume* pv;


    
    //-- Get all disk /crystal informations here
    //   crystal longitudinal size are defined w.r.t the center of the hexagon (like a radius)
    //   crystal/disk length are defined the usual way (each wrapper adds 2*wrapping_size!)
        
    DiskCalorimeter const & cal       = *(GeomHandle<DiskCalorimeter>());

    //calorimeter mother enveloppe
    G4double mother_inRadius          = cal.caloGeomInfo().enveloppeInRadius();
    G4double mother_outRadius         = cal.caloGeomInfo().enveloppeOutRadius();
    G4double mother_z0                = cal.caloGeomInfo().enveloppeZ0();
    G4double mother_z1                = cal.caloGeomInfo().enveloppeZ1();
    
    //crystal properties
    G4int    nRO                      = cal.caloGeomInfo().nROPerCrystal();
    G4double ROHalfThickness          = cal.caloGeomInfo().roHalfThickness();
    G4double ROHalfTrans              = cal.caloGeomInfo().roHalfTrans();

    G4int    crystalnEdges            = cal.caloGeomInfo().crystalNedges();
    G4double crystalPolysize          = cal.caloGeomInfo().crystalHalfTrans();
    G4double crystalDepth             = 2.0*cal.caloGeomInfo().crystalHalfLength();

    G4double wrapThickness            = cal.caloGeomInfo().wrapperThickness();
    G4double wrapPolysize             = crystalPolysize + wrapThickness;
    G4double wrapDepth                = crystalDepth + 2.0*ROHalfThickness + 2.0*wrapThickness; 
    
    //calorimeter calibration system
    G4int nPipes                      = cal.caloGeomInfo().nPipes();      
    G4double pipeRadius               = cal.caloGeomInfo().pipeRadius();
    G4double pipeThickness            = cal.caloGeomInfo().pipeThickness();
    std::vector<double> pipeTorRadius = cal.caloGeomInfo().pipeTorRadius();

    //disk properties
    const unsigned int nDisks         = cal.nDisk();
    G4double caseThickness            = cal.caloGeomInfo().caseThickness();
    G4double caseDepth                = wrapDepth + 2.0*caseThickness;
    G4double diskDepth                = caseDepth + 2.0*pipeRadius;

    
    // Readout positions
    std::vector<double> XposRO, YposRO;

    double R0disp = 0.5*crystalPolysize;
    if (nRO==1) {XposRO.push_back(0);
                 YposRO.push_back(0);}
    if (nRO==2) {XposRO.push_back(0);XposRO.push_back(0); 
                 YposRO.push_back(-R0disp);YposRO.push_back(R0disp);}
    if (nRO==4) {XposRO.push_back(-R0disp);XposRO.push_back(-R0disp);XposRO.push_back(R0disp);XposRO.push_back(R0disp);
                 YposRO.push_back(-R0disp);YposRO.push_back(R0disp);YposRO.push_back(-R0disp);YposRO.push_back(R0disp);}






    // crystal z position
    G4double ZPoscrystal  = wrapThickness;
    G4double ZPosR0       = wrapThickness+crystalDepth+ROHalfThickness;


    // disk position
    G4double ZposPipe     = -diskDepth/2.0+pipeRadius;                                            
    G4double ZposCase     = pipeRadius;                                            









    //--------------------------------------
    // Building blocks for a crystal
    //

    //
    // define required solids
    double offsetAngle = (crystalnEdges==4) ? CLHEP::pi/4.0 : 0;  //need to add an offset to the phi angle to _squares_ to have them rotated properly

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

    G4Box *crystalRO = new G4Box("CrystalRO",ROHalfTrans,ROHalfTrans,ROHalfThickness);



    //
    // define required logical volumes

    G4LogicalVolume *CrystalLog  = new G4LogicalVolume(crystal, crysMaterial, "CrystalLog");
    G4LogicalVolume *ROLog       = new G4LogicalVolume(crystalRO, readMaterial, "CrystalROLog" );    
    
    G4VisAttributes* crys_visAtt = new G4VisAttributes(isCrystalVisible, G4Color::Green());
    crys_visAtt->SetForceSolid(isCrystalSolid);
    crys_visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);

    CrystalLog->SetVisAttributes(crys_visAtt);    
    ROLog->SetVisAttributes(G4VisAttributes::Invisible);
   
    //-- Sensitive detector
    G4VSensitiveDetector* ccSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloCrystal());
    G4VSensitiveDetector* crSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloReadout());
    
    CrystalLog->SetSensitiveDetector(ccSD);
    ROLog->SetSensitiveDetector(crSD);







  

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
	double zhl         = static_cast<G4Tubs*>(calorimeterInfo.solid)->GetZHalfLength();
	CLHEP::Hep3Vector const & CalorimeterOffsetInMu2e = calorimeterInfo.centerInMu2e();
	double CalorimeterOffsetInMu2eZ = CalorimeterOffsetInMu2e[CLHEP::Hep3Vector::Z];
	cout << __func__ << " Calorimeter mother center in Mu2e   : " << CalorimeterOffsetInMu2e << endl;
	cout << __func__ << " Calorimeter mother Z extent in Mu2e : " << CalorimeterOffsetInMu2eZ - zhl << ", " << CalorimeterOffsetInMu2eZ + zhl << endl;
    }












    //--------------------------------------
    // Construct disks: diskCase contains the crystals. DiskBox contains the calibration pipes and diskCase
    //
    VolumeInfo diskBoxInfo[nDisks];
    VolumeInfo diskCaseInfo[nDisks];
    int nTotCrystal(0);




    for (unsigned int idisk=0;idisk<nDisks;++idisk)
    {

	ostringstream discname0;      discname0<<"DiskCalorimeter_" <<idisk;
	ostringstream discname1;      discname1<<"DiskCase_" <<idisk;

	double radiusIn       = cal.disk(idisk).innerRadius()-caseThickness;
	double radiusOut      = cal.disk(idisk).outerRadius()+caseThickness;
	double diskpar0[5]    = {radiusIn,radiusOut, diskDepth/2.0, 0., CLHEP::twopi};
	double diskpar1[5]    = {radiusIn,radiusOut, caseDepth/2.0, 0., CLHEP::twopi};

	G4ThreeVector posDisk = cal.disk(idisk).origin() - posCaloMother;


	diskBoxInfo[idisk] =  nestTubs(discname0.str(),
                        	       diskpar0,fillMaterial,&cal.disk(idisk).rotation(),posDisk,
                        	       calorimeterInfo,
                        	       idisk,
                        	       isDiskCaseVisible,G4Colour::Green(),isDiskCaseSolid,forceAuxEdgeVisible,
                        	       true,doSurfaceCheck );

	diskCaseInfo[idisk] = nestTubs(discname1.str(),
                        	       diskpar1,diskMaterial,0,G4ThreeVector(0.0,0.0,ZposCase),
                        	       diskBoxInfo[idisk],
                        	       nTotCrystal,
                        	       isDiskCaseVisible,G4Colour::Red(),isDiskCaseSolid,forceAuxEdgeVisible,
                        	       true,doSurfaceCheck );
			      
	if ( verbosityLevel > 0) 
	{
	    cout << __func__ << " CalorimeterDisk center in Mu2e    : " << diskBoxInfo[idisk].centerInMu2e() << endl;
	    cout << __func__ << " CalorimeterDisk Z extent in Mu2e  : " << diskBoxInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] - diskDepth/2.0 << ", " << diskBoxInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] + diskDepth/2.0 << endl;
	    cout << __func__ << " CalorimeterCase center in Mu2e    : " << diskCaseInfo[idisk].centerInMu2e() << endl;
	    cout << __func__ << " CalorimeterCase Z extent in Mu2e  : " << diskCaseInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] - caseDepth/2.0 << ", " << diskCaseInfo[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] + caseDepth/2.0 << endl;
	}


	// fill it with crystals
	for(int ic=0; ic <cal.disk(idisk).nCrystals(); ++ic)
	{

	      G4int id       = nTotCrystal+ic;
	      G4int roidBase = cal.ROBaseByCrystal(id);

              //position of shell in the disk
              CLHEP::Hep3Vector crystalPosition = cal.disk(idisk).crystal(ic).localPosition();
              double x = crystalPosition.x();
              double y = crystalPosition.y();
              double z = -wrapDepth/2.0; 	      


	      G4LogicalVolume *thisWrapLog = new G4LogicalVolume(crystalWrap, wrapMaterial, "WrapLog");
	      thisWrapLog->SetVisAttributes(G4VisAttributes::Invisible);
	      pv = new G4PVPlacement(0,G4ThreeVector(x,y,z),thisWrapLog,"CrysWrapPV",diskCaseInfo[idisk].logical,false,id,false);
              doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

              pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,ZPoscrystal),CrystalLog,"CrysPV",thisWrapLog,false,id,false);
              doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

              // add the readout
	      for (unsigned int iro=0;iro < XposRO.size();++iro)
	      {
		 pv = new G4PVPlacement(0,G4ThreeVector(XposRO[iro],YposRO[iro],ZPosR0),ROLog,"CrysROPV_0",thisWrapLog,true,roidBase+iro,false);
		 doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
	      }

	   }
   
           nTotCrystal += cal.disk(idisk).nCrystals();
     }




     //--------------------------------------
     // Construct Calibration system
     //
     if (nPipes>0) 
     {
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
