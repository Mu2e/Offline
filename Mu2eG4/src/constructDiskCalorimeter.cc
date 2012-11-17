//
// Free function to create the calorimeter.
//
//
// Original author Ivan Logashenko
//
// Notes
// 1) The argument zOff is the zlocation of the center of the mother volume,
//    as mesaured in the mu2e coordinate system.
// 2) The vanes are placed directly into DS3.  We did not make a mother volume for them.
//
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
#include "Mu2eG4/inc/constructDiskCalorimeter.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"
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
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"




using namespace std;

namespace mu2e {

  void constructDiskCalorimeter( VolumeInfo const &  mother,
                                 double              zOffset,
                                 SimpleConfig const& config ){



    //-- Read parameters from config file
    int const verbosityLevel        = config.getInt("calorimeter.verbosityLevel",0);
    bool const isDiskBoxVisible     = config.getBool("calorimeter.vaneBoxVisible",true);
    bool const isDiskBoxSolid       = config.getBool("calorimeter.vaneBoxSolid",true);
    bool const isCrystalVisible     = config.getBool("calorimeter.crystalVisible",false);
    bool const isCrystalSolid       = config.getBool("calorimeter.crystalSolid",true);
    bool const forceAuxEdgeVisible  = config.getBool("g4.forceAuxEdgeVisible",false);
    G4bool const doSurfaceCheck     = config.getBool("g4.doSurfaceCheck",false);
    bool const placePV              = true;



    //-- A helper class for parsing the config file.
    MaterialFinder materialFinder(config);

    G4Material* diskMaterial = materialFinder.get("calorimeter.calorimeterDiskMaterial");
    G4Material* fillMaterial = materialFinder.get("calorimeter.calorimeterFillMaterial");
    G4Material* crysMaterial = materialFinder.get("calorimeter.crystalMaterial");
    G4Material* wrapMaterial = materialFinder.get("calorimeter.crystalWrapper");
    G4Material* readMaterial = materialFinder.get("calorimeter.crystalReadoutMaterial");
    
    
    
    //-- Get all disk /crystal informations here
    //   crystal longitudinal size are defined w.r.t the center of the hexagon (like a radius)
    //   crystal/disk length are defined the usual way (each wrapper adds 2*wrapping_size!)
        
    DiskCalorimeter const & cal    = *(GeomHandle<DiskCalorimeter>());

    G4int    nRO                   = cal.nROPerCrystal();
    G4double ROHalfThickness       = cal.roHalfThickness();
    G4double ROHalfTrans           = cal.roHalfSize();

    G4double crystalHexsize        = cal.crysHalfTrans();
    G4double crystalDepth          = cal.crysDepth();

    G4double wrapThickness         = cal.wrapperThickness();
    G4double wrapHexsize           = crystalHexsize + wrapThickness;
    G4double wrapDepth             = crystalDepth + 2.0*ROHalfThickness + 2.0*wrapThickness; 
    
    G4double shellThickness        = cal.shellThickness();
    G4double shellHexsize          = wrapHexsize + shellThickness;
    G4double shellDepth            = wrapDepth; //+ 2.0*shellThickness; for a crystal shell surrounding the z faces

    
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
    G4Box *crystalRO              = new G4Box("CrystalRO",ROHalfTrans,ROHalfTrans,ROHalfThickness);






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














    G4ThreeVector pcalo = cal.getOrigin();

    //-- Construct disks: diskInner contains the crystals, and is inside diskBox

    const unsigned int nDisks = cal.nDisks();
    VolumeInfo diskBoxInfo[nDisks];
    VolumeInfo diskInnerInfo[nDisks];

    //counter of crystals in disks
    G4int crystalIdOffset(0);
    
    for (unsigned int idisk=0;idisk<nDisks;++idisk){

        G4ThreeVector pdisk = cal.getDisk(idisk).getOriginLocal();
double hack = 600;	
	G4ThreeVector pos = G4ThreeVector(pdisk.x(),pdisk.y(), pcalo.z()+pdisk.z()+zOffset-hack);

	ostringstream discname1;      discname1<<"DiskCalorimeter_" <<idisk;
	ostringstream discname2;      discname2<<"DiskInner_" <<idisk;

	double radiusIn   = cal.getDisk(idisk).innerRadius();
	double radiusOut  = cal.getDisk(idisk).outerRadius();
	double diskDepth  = shellDepth + 2.0*cal.getDisk(idisk).thickness();

	double diskpar1[5] = {radiusIn-cal.getDisk(idisk).thickness(),radiusOut+cal.getDisk(idisk).thickness(), diskDepth/2.0, 0, 2*pi};
	double diskpar2[5] = {radiusIn,radiusOut, shellDepth/2.0, 0, 2*pi};


	diskBoxInfo[idisk] = nestTubs(discname1.str(),
                              diskpar1,
                              diskMaterial,
                              &cal.getDisk(idisk).getRotation(),
                              pos,
                              mother,
                              idisk,
                              isDiskBoxVisible,
                              G4Colour::Yellow(),
                              isDiskBoxSolid,
                              forceAuxEdgeVisible,
                              placePV,
                              doSurfaceCheck );

	diskInnerInfo[idisk] = nestTubs(discname2.str(),
                              diskpar2,
                              fillMaterial,
                              0,
                              G4ThreeVector(0.0,0.0,0.0),
                              diskBoxInfo[idisk],
                              100*idisk,
                              isDiskBoxVisible,
                              G4Colour::Magenta(),
                              isDiskBoxSolid,
                              forceAuxEdgeVisible,
                              placePV,
                              doSurfaceCheck );
			      

	if ( verbosityLevel > 0) {
	    double zhl         = diskDepth/2.0;
	    CLHEP::Hep3Vector const & CalorimeterDiskOffsetInMu2e = diskBoxInfo[idisk].centerInMu2e();
	    double CalorimeterDiskOffsetInMu2eZ = CalorimeterDiskOffsetInMu2e[CLHEP::Hep3Vector::Z];
	    cout << __func__ << " CalorimeterDisk center in Mu2e      : " << CalorimeterDiskOffsetInMu2e << endl;
	    cout << __func__ << " CalorimeterDisk Z extent in Mu2e    : " <<
            CalorimeterDiskOffsetInMu2eZ - zhl << ", " << CalorimeterDiskOffsetInMu2eZ + zhl << endl;
	}



	// -- then fill the inner disk with crystals
	G4int nCrystalInThisDisk = cal.getDisk(idisk).nCrystals();			
	for(int ic=0; ic <nCrystalInThisDisk; ++ic){

	      G4int id       = crystalIdOffset+ic;
	      G4int roidBase = cal.getROBaseByCrystal(id);


	      // Have to define a shell / wrapper logical volume for each crystal 
	      // to get correct index in CrystalCaloSD
	      G4LogicalVolume *thisShellLog(0);
	      if (shellThickness > 0.001){
		thisShellLog = new G4LogicalVolume(crystalShell, fillMaterial, "ShellLog");
		thisShellLog->SetVisAttributes(G4VisAttributes::Invisible);
	      }

	      G4LogicalVolume *thisWrapLog = new G4LogicalVolume(crystalWrap, wrapMaterial, "WrapLog");
	      thisWrapLog->SetVisAttributes(G4VisAttributes::Invisible);

              //position of shell in the disk
	      //contrary to rectangles, z position of hexagon is their base, not their center in Geant 4!!	      
              CLHEP::Hep3Vector crystalPosition = cal.getDisk(idisk).getCrystal(ic).position();
              double x = crystalPosition.x();
              double y = crystalPosition.y();
              double z = -shellDepth/2.0; 	      

              // place a shell only if it has non-zero thickness, or place the wrapper directly
              if (shellThickness > 0.001) {
	           new G4PVPlacement(0,G4ThreeVector(x,y,z),thisShellLog,"CrysShellPV",diskInnerInfo[idisk].logical,0,id,doSurfaceCheck);   
                   new G4PVPlacement(0,G4ThreeVector(0.0,0.0,wrapZPos),thisWrapLog,"CrysWrapPV",thisShellLog,0,id,doSurfaceCheck);
              } else {
	           new G4PVPlacement(0,G4ThreeVector(x,y,z),thisWrapLog,"CrysWrapPV",diskInnerInfo[idisk].logical,0,id,doSurfaceCheck);   	      
	      }

	      // -- place crystal inside warp
              new G4PVPlacement(0,G4ThreeVector(0.0,0.0,crystalZPos),CrystalLog,"CrysPV",thisWrapLog,0,id,doSurfaceCheck);


              // -- add the readout
	      if (nRO==1) 
		new G4PVPlacement(0,G4ThreeVector(0,0,ROZPos),ROLog,"CrysROPV_0",thisWrapLog,0,roidBase,doSurfaceCheck);

	      if (nRO==2) { 
		new G4PVPlacement(0,G4ThreeVector(0,-0.5*crystalHexsize,ROZPos),ROLog,"CrysROPV_0",thisWrapLog,0,roidBase,doSurfaceCheck);
        	new G4PVPlacement(0,G4ThreeVector(0, 0.5*crystalHexsize,ROZPos),ROLog,"CrysROPV_1",thisWrapLog,0,roidBase+1,doSurfaceCheck);
              }

	      if (nRO==4) { 
		new G4PVPlacement(0,G4ThreeVector(-0.5*crystalHexsize,-0.5*crystalHexsize,ROZPos),ROLog,"CrysROPV_0",thisWrapLog,0,roidBase,doSurfaceCheck);
        	new G4PVPlacement(0,G4ThreeVector(-0.5*crystalHexsize, 0.5*crystalHexsize,ROZPos),ROLog,"CrysROPV_1",thisWrapLog,0,roidBase+1,doSurfaceCheck);
        	new G4PVPlacement(0,G4ThreeVector( 0.5*crystalHexsize,-0.5*crystalHexsize,ROZPos),ROLog,"CrysR0PV_2",thisWrapLog,0,roidBase+2,doSurfaceCheck);
        	new G4PVPlacement(0,G4ThreeVector( 0.5*crystalHexsize, 0.5*crystalHexsize,ROZPos),ROLog,"CrysROPV_3",thisWrapLog,0,roidBase+3,doSurfaceCheck);
              }


	   }//end loop over crystals
           
	   crystalIdOffset +=nCrystalInThisDisk;

     }//end loop over disks


  

  }//end of disk calo construction

} // end namespace mu2e
