//
// Free function to create the Vane calorimeter.
//
// $Id: constructVaneCalorimeter.cc,v 1.14 2014/08/01 23:14:36 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 23:14:36 $
//
// Original author Ivan Logashenko
// Modified by Bertrand Echenard
//
// Notes
//
//  1. a crystal has readouts at the back, both are surrounded by the wrapping, and the wrapping by a shell
//  2. by default, the wrapping surrounds the front/back face of the crystal+ro, the shell does not (shell is a casing)
//  3. The vanes are placed directly into DS3.  We did not make a mother volume for them.
//  4. The argument zOff is the zlocation of the center of the mother volume, as mesaured in the mu2e coordinate system.
//
//  5) Modified version  builds the calorimeter by making a physical mother volume.

#include <iostream>

// Mu2e includes.
#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"
#include "Mu2eG4/inc/constructVaneCalorimeter.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/Vane.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "Mu2eG4/inc/CaloCrystalSD.hh"
#include "Mu2eG4/inc/CaloReadoutSD.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"

// G4 includes
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4Tubs.hh"

using namespace std;

namespace mu2e {

  VolumeInfo constructVaneCalorimeter( VolumeInfo const &  mother,SimpleConfig const& config )
  {


    //-- Read parameters from config file
    //   Control of graphics for debugging the geometry.
    //   Only instantiate sectors to be drawn.
    int const verbosityLevel             = config.getInt("calorimeter.verbosityLevel",0);


    bool const isCalorimeterVisible      = config.getBool("calorimeter.calorimeterVisible",false);
    bool const isCalorimeterSolid        = config.getBool("calorimeter.calorimeterSolid",false);
    bool const isVaneBoxVisible          = config.getBool("calorimeter.boxVisible",true);
    bool const isVaneBoxSolid            = config.getBool("calorimeter.boxSolid",true);
    bool const isAbsorberBoxVisible      = config.getBool("calorimeter.absorberBoxVisible",true);
    bool const isAbsorberBoxSolid        = config.getBool("calorimeter.absorberBoxSolid",true);
    bool const isCrystalVisible          = config.getBool("calorimeter.crystalVisible",false);
    bool const isCrystalSolid            = config.getBool("calorimeter.crystalSolid",true);

    bool const forceAuxEdgeVisible       = config.getBool("g4.forceAuxEdgeVisible",false);
    G4bool const doSurfaceCheck          = config.getBool("g4.doSurfaceCheck",false);
    bool const placePV                   = true;


    //calorimeter mother neveloppe
    double mother_inRadius  = config.getDouble("calorimeter.caloMotherInRadius",0); 
    double mother_outRadius = config.getDouble("calorimeter.caloMotherOutRadius",765); 
    double mother_z0        = config.getDouble("calorimeter.caloMotherZ0",11740); 
    double mother_z1        = config.getDouble("calorimeter.caloMotherZ1",13910); 


    //-- A helper class for parsing the config file.
    MaterialFinder materialFinder(config);
    G4Material* fillMaterial             = materialFinder.get("calorimeter.calorimeterFillMaterial");
    G4Material* crysMaterial             = materialFinder.get("calorimeter.crystalMaterial");
    G4Material* wrapMaterial             = materialFinder.get("calorimeter.crystalWrapper");
    G4Material* readMaterial             = materialFinder.get("calorimeter.crystalReadoutMaterial");
    G4Material* shieldMaterial           = materialFinder.get("calorimeter.shieldMaterial");
    G4Material* neutronAbsorberMaterial  = materialFinder.get("calorimeter.neutronAbsorberMaterial");


    //-- parameters for beam scintillator trigger for TestBeam configuration
    int        hasScint                  = config.getInt("calorimeter.hasScint", 0);
    G4Material*scintMaterial             = materialFinder.get("calorimeter.scintMaterial", "G4_PLASTIC_SC_VINYLTOLUENE");
    double     scintDx                   = config.getDouble("calorimeter.scintDz", 20.);
    double     scintDy                   = config.getDouble("calorimeter.scintDz", 20.);
    double     scintDz                   = config.getDouble("calorimeter.scintDz", 10.);
    std::vector<double> scintPos;
    config.getVectorDouble("calorimeter.scintPos", scintPos, 3);


    //-- Get calorimeter handle
    VaneCalorimeter const & cal = *(GeomHandle<VaneCalorimeter>());






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
      cout << __func__ << " Calorimeter mother Z extent in Mu2e    : " <<CalorimeterOffsetInMu2eZ - zhl << ", " << CalorimeterOffsetInMu2eZ + zhl << endl;
    }

    G4int nRO                   = cal.caloGeomInfo().nROPerCrystal();
    G4double crystalSize        = cal.caloGeomInfo().crystalHalfTrans();
    G4double crystalLength      = cal.caloGeomInfo().crystalHalfLength();
    G4double wrapSize           = crystalSize   + cal.caloGeomInfo().wrapperThickness();
    G4double wrapLength         = crystalLength + cal.caloGeomInfo().wrapperThickness() + cal.caloGeomInfo().roHalfThickness();
    G4double shellSize          = wrapSize      + cal.caloGeomInfo().shellThickness();
    G4double shellLength        = wrapLength;


    //-- Create solids for one crystal
    G4Box *crystalShell = new G4Box("CrystalShell",shellLength,shellSize,shellSize);
    G4Box *crystalWrap  = new G4Box("CrystalWrap",wrapLength,wrapSize,wrapSize);
    G4Box *crystal      = new G4Box("Crystal",crystalLength,crystalSize,crystalSize);
    G4Box *crystalRO    = new G4Box("CrystalRO",cal.caloGeomInfo().roHalfThickness(),cal.caloGeomInfo().roHalfTrans(),cal.caloGeomInfo().roHalfTrans() );


    //-- Definition of a few logical volumes    
    //
    // Geant4 indexing is such that the crystal /readout can be defined once, 
    // but the shell and wrapper logical volumes must be defined every time a crystal is placed  
    // to get correct index in CaloCrystalSD class

    G4LogicalVolume *CrystalLog  = new G4LogicalVolume(crystal, crysMaterial, "CrystalLog");
    G4LogicalVolume *ROLog       = new G4LogicalVolume(crystalRO, readMaterial, "CrystalROLog" );    

    if(!isCrystalVisible) 
    {
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





    //-- Build vanes and absorbers
    const int nvane = cal.nVane();
    const double shieldHalfThickness           = config.getDouble("calorimeter.shieldHalfThickness");
    const double neutronAbsorberHalfThickness  = config.getDouble("calorimeter.neutronAbsorberHalfThickness");
    const double absorberHalfThickness         = shieldHalfThickness + neutronAbsorberHalfThickness;


    VolumeInfo vaneOutInfo[nvane];
    VolumeInfo vaneInInfo[nvane];
    VolumeInfo shieldInfo[nvane];
    VolumeInfo neutronAbsorberInfo[nvane];
    
    //-- plastic scintillator finger for test beam
    VolumeInfo scintInfo;

    G4VPhysicalVolume* pv;


    for( int iv=0; iv<nvane; ++iv ) 
    {

      ostringstream nameOutVane;           nameOutVane          << "CalorimeterOutVane_"            << iv;
      ostringstream nameInVane;            nameInVane           << "CalorimeterInVane_"             << iv;
      ostringstream nameShield;            nameShield           << "CalorimeterShield_"             << iv;
      ostringstream nameNeutronAbsorber;   nameNeutronAbsorber  << "CalorimeterNeutronAbsorber_"    << iv;
      
      const CLHEP::Hep3Vector & sizeOut = cal.vane(iv).size();
      const CLHEP::Hep3Vector & sizeIn  = cal.vane(iv).size() - CLHEP::Hep3Vector(cal.caloGeomInfo().caseThickness(),cal.caloGeomInfo().caseThickness(),cal.caloGeomInfo().caseThickness());
      
      double dimOutVane[3]             = {sizeOut.x(), sizeOut.y(), sizeOut.z()};
      double dimInVane[3]              = {sizeIn.x(),  sizeIn.y(),  sizeIn.z() - absorberHalfThickness};      
      double dimShield[3]              = {sizeIn.x(),  sizeIn.y(),  shieldHalfThickness};
      double dimNeutronAbsorber[3]     = {sizeIn.x(),  sizeIn.y(),  neutronAbsorberHalfThickness};
 
      G4ThreeVector posVane            = cal.origin() + cal.vane(iv).originLocal() - posCaloMother;
      G4ThreeVector posInVane          = G4ThreeVector(0., 0. , absorberHalfThickness);
      G4ThreeVector posShield          = G4ThreeVector(0., 0., -sizeIn.z() + shieldHalfThickness);
      G4ThreeVector posNeutronAbsorber = G4ThreeVector(0., 0. ,-sizeIn.z() + 2*shieldHalfThickness + neutronAbsorberHalfThickness);


      vaneOutInfo[iv]     = nestBox(nameOutVane.str(),
				 dimOutVane,
				 fillMaterial,
				 &cal.vane(iv).rotation(),
				 posVane,
				 calorimeterInfo,
				 iv,
				 isVaneBoxVisible,
				 G4Colour::Yellow(),
				 isVaneBoxSolid,
				 forceAuxEdgeVisible,
				 placePV,
				 doSurfaceCheck );


      vaneInInfo[iv]     = nestBox(nameInVane.str(),
				 dimInVane,
				 fillMaterial,
				 0,
				 posInVane,
				 vaneOutInfo[iv],
				 iv,
				 isVaneBoxVisible,
				 G4Colour::Yellow(),
				 isVaneBoxSolid,
				 forceAuxEdgeVisible,
				 placePV,
				 doSurfaceCheck );

      if (verbosityLevel) cout << "Calorimeter Vane position: "<<posVane<<endl;

      
      if( shieldHalfThickness > 0.0)
      {
	shieldInfo[iv]    = nestBox(nameShield.str(),
				    dimShield,
				    shieldMaterial,
				    0,
				    posShield,
				    vaneOutInfo[iv] ,
				    iv,
				    isAbsorberBoxVisible,
				    G4Colour::Blue(),
				    isAbsorberBoxSolid,
				    forceAuxEdgeVisible,
				    placePV,
				    doSurfaceCheck );
      }

      if( neutronAbsorberHalfThickness > 0)
      {
	neutronAbsorberInfo[iv]    = nestBox(nameNeutronAbsorber.str(),
					     dimNeutronAbsorber,
					     neutronAbsorberMaterial,
					     0,
					     posNeutronAbsorber,
					     vaneOutInfo[iv] ,
					     iv,
					     isAbsorberBoxVisible,
					     G4Colour::Cyan(),
					     isAbsorberBoxSolid,
					     forceAuxEdgeVisible,
					     placePV,
					     doSurfaceCheck );
      }
      
      if ( verbosityLevel > 0) 
      {
	double xhl  = static_cast<G4Box*>(vaneOutInfo[iv].solid)->GetXHalfLength();
	double yhl  = static_cast<G4Box*>(vaneOutInfo[iv].solid)->GetYHalfLength();
	double zhl  = static_cast<G4Box*>(vaneOutInfo[iv].solid)->GetZHalfLength();
	cout << __func__ << " center in Mu2e    : " << vaneOutInfo[iv].centerInMu2e() << endl;
	cout << __func__ << " X extent in Mu2e  : " <<vaneOutInfo[iv].centerInMu2e().x() - xhl << ", " << vaneOutInfo[iv].centerInMu2e().x() + xhl << endl;
	cout << __func__ << " Y extent in Mu2e  : " <<vaneOutInfo[iv].centerInMu2e().y() - yhl << ", " << vaneOutInfo[iv].centerInMu2e().y() + yhl << endl;
	cout << __func__ << " Z extent in Mu2e  : " <<vaneOutInfo[iv].centerInMu2e().z() - zhl << ", " << vaneOutInfo[iv].centerInMu2e().z() + zhl << endl;



      }




      //-- place crystals inside vanes
      G4int ncrys                 = cal.vane(iv).nCrystals();
      for( int ic=0; ic<ncrys; ++ic ) 
      {

	// IDs
	G4int id       = iv*ncrys + ic;       // Crystal ID
	G4int roidBase = cal.ROBaseByCrystal(id);


	// Have to define a shell / wrapper logical volume for each crystal 
	// to get correct index in CrystalCaloSD
	G4LogicalVolume *thisShellLog(0);
	if (cal.caloGeomInfo().shellThickness() > 0.001)
	{
	  thisShellLog = new G4LogicalVolume(crystalShell, fillMaterial, "ShellLog");
	  thisShellLog->SetVisAttributes(G4VisAttributes::Invisible);
	}

	G4LogicalVolume *thisWrapLog = new G4LogicalVolume(crystalWrap, wrapMaterial, "WrapLog");
	thisWrapLog->SetVisAttributes(G4VisAttributes::Invisible);



	// Position - first run along Z, then along Y, both times in positive direction
	//this is the position of the wrapper, so the x must be zero, not crystalPosition.x()
	CLHEP::Hep3Vector crystalPosition = cal.vane(iv).crystal(ic).localPosition();
	double x = 0;
	double y = crystalPosition.y();
	double z = crystalPosition.z();
	     

	// place a shell only if it has non-zero thickness, or place the wrapper directly
	if (cal.caloGeomInfo().shellThickness()  > 0.001) 
	{
	  pv = new G4PVPlacement(0,G4ThreeVector(x,y,z),thisShellLog,"CrysShellPV",
                            vaneInInfo[iv].logical,false,id,false);
          doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

	  pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,0.0),thisWrapLog,"CrysWrapPV",
                            thisShellLog,false,id,false);
          doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

	} else {
	  pv = new G4PVPlacement(0,G4ThreeVector(x,y,z),thisWrapLog,"CrysWrapPV",
                            vaneInInfo[iv].logical,false,id,false);   	      
          doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

	}

//swapped sign from all these

	// -- place crystal inside warp
	pv = new G4PVPlacement(0,G4ThreeVector(-cal.caloGeomInfo().roHalfThickness(),0.0,0.0),CrystalLog,"CrysPV",
                          thisWrapLog,false,id,false);
        doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

	// -- add the readout
	if (nRO==1) 
	  pv = new G4PVPlacement(0,G4ThreeVector(cal.caloGeomInfo().crystalHalfLength(),0,0),ROLog,"CrysROPV_0",
                            thisWrapLog,false,roidBase,false);
          doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

	if (nRO==2) 
	{ 
	  pv = new G4PVPlacement(0,G4ThreeVector(cal.caloGeomInfo().crystalHalfLength(),-0.5*cal.caloGeomInfo().crystalHalfTrans(),0.0),ROLog,"CrysROPV_0",
                            thisWrapLog,false,roidBase,false);
          doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
	  pv = new G4PVPlacement(0,G4ThreeVector(cal.caloGeomInfo().crystalHalfLength(), 0.5*cal.caloGeomInfo().crystalHalfTrans(),0.0),ROLog,"CrysROPV_1",
                            thisWrapLog,false,roidBase+1,false);
          doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
	}

	if (nRO==4) 
	{ 
	  G4double cHS = -0.5*cal.caloGeomInfo().crystalHalfTrans();
	  pv = new G4PVPlacement(0,G4ThreeVector(cal.caloGeomInfo().crystalHalfLength(),-cHS,-cHS),ROLog,"CrysROPV_0"
                            ,thisWrapLog,false,roidBase,false);
          doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
	  pv = new G4PVPlacement(0,G4ThreeVector(cal.caloGeomInfo().crystalHalfLength(),-cHS,cHS), ROLog,"CrysROPV_1",
                            thisWrapLog,false,roidBase+1,false);
          doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
	  pv = new G4PVPlacement(0,G4ThreeVector(cal.caloGeomInfo().crystalHalfLength(), cHS,-cHS),ROLog,"CrysR0PV_2",
                            thisWrapLog,false,roidBase+2,false);
          doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
          pv = new G4PVPlacement(0,G4ThreeVector(cal.caloGeomInfo().crystalHalfLength(), cHS,cHS), ROLog,"CrysROPV_3",
                            thisWrapLog,false,roidBase+3,false);
          doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
	}


      }//end crystal loop



    }//end loop over vanes

    if (hasScint == 1){
      double posx = scintPos[0] - posCaloMother.x();
      double posy = scintPos[1] - posCaloMother.y();
      double posz = scintPos[2] - posCaloMother.z();
      G4ThreeVector pos(posx, posy, posz);
      double dimScint[3] = {scintDx, scintDy, scintDz};

      scintInfo =  nestBox("ScintillatorFinger",
			   dimScint,
			   scintMaterial,
			   0,
			   pos,
			   calorimeterInfo,
			   0,
			   true,
			   G4Colour::Red(),
			   isVaneBoxSolid,
			   forceAuxEdgeVisible,
			   placePV,
			   doSurfaceCheck );
    }

    return calorimeterInfo;

  }


} // end namespace mu2e
