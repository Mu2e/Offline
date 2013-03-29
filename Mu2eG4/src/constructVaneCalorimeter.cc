//
// Free function to create the calorimeter.
//
// $Id: constructVaneCalorimeter.cc,v 1.5 2013/03/29 04:35:17 gandr Exp $
// $Author: gandr $
// $Date: 2013/03/29 04:35:17 $
//
// Original author Ivan Logashenko
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

  VolumeInfo constructVaneCalorimeter( VolumeInfo const &  mother,
                                       double zOffset,
                                       SimpleConfig const& config ){

    // Access to the G4HelperService.
    G4Helper    & _helper = *(art::ServiceHandle<G4Helper>());


    //-- Read parameters from config file
    //   Control of graphics for debugging the geometry.
    //   Only instantiate sectors to be drawn.
    int const verbosityLevel             = config.getInt("calorimeter.verbosityLevel",0);


    bool const isCalorimeterVisible      = config.getBool("calorimeter.calorimeterVisible",true);
    bool const isCalorimeterSolid        = config.getBool("calorimeter.calorimeterSolid",true);
    bool const isVaneBoxVisible          = config.getBool("calorimeter.vaneBoxVisible",true);
    bool const isVaneBoxSolid            = config.getBool("calorimeter.vaneBoxSolid",true);
    bool const isAbsorberBoxVisible      = config.getBool("calorimeter.absorberBoxVisible",true);
    bool const isAbsorberBoxSolid        = config.getBool("calorimeter.absorberBoxSolid",true);
    bool const isCrystalVisible          = config.getBool("calorimeter.crystalVisible",false);
    bool const isCrystalSolid            = config.getBool("calorimeter.crystalSolid",true);
    bool const isShieldVisible           = config.getBool("calorimeter.shieldVisible",false);
    bool const isShieldSolid             = config.getBool("calorimeter.shieldSolid",true);
    bool const isNeutronAbsorberVisible  = config.getBool("calorimeter.neutronAbsorberVisible",false);
    bool const isNeutronAbsorberSolid    = config.getBool("calorimeter.neutronAbsorberSolid",true);

    bool const forceAuxEdgeVisible       = config.getBool("g4.forceAuxEdgeVisible",false);
    G4bool const doSurfaceCheck          = config.getBool("g4.doSurfaceCheck",false);
    bool const placePV                   = true;




    //-- A helper class for parsing the config file.
    MaterialFinder materialFinder(config);
    G4Material* fillMaterial             = materialFinder.get("calorimeter.calorimeterFillMaterial");
    G4Material* crysMaterial             = materialFinder.get("calorimeter.crystalMaterial");
    G4Material* wrapMaterial             = materialFinder.get("calorimeter.crystalWrapper");
    G4Material* readMaterial             = materialFinder.get("calorimeter.crystalReadoutMaterial");
    G4Material* shieldMaterial           = materialFinder.get("calorimeter.shieldMaterial");
    G4Material* neutronAbsorberMaterial  = materialFinder.get("calorimeter.neutronAbsorberMaterial");

    //-- Get calorimeter handle
    VaneCalorimeter const & cal = *(GeomHandle<VaneCalorimeter>());
    double rIn                           = cal.envelopeRmin();
    double rOut                          = cal.envelopeRmax();
    double halfLength                    = cal.envelopeHalfLength();

    //  Make the mother volume for the calorimeter.
    G4ThreeVector pcalo = cal.origin();
    G4ThreeVector pos   = G4ThreeVector(0, 0, pcalo.z()+zOffset);
    verbosityLevel > 0 &&
      cout << "Calorimeter  position Origin local: ("
	   << pos.x() << "," << pos.y() << "," << pos.z()
	   << ")" << endl;
	      

    TubsParams caloParams( rIn, rOut, halfLength, 0., CLHEP::twopi);

    VolumeInfo calorimeterInfo = nestTubs( "CalorimeterMother",
					   caloParams,
					   fillMaterial,
					   0,
					   pos,
					   mother,
					   0,
					   isCalorimeterVisible,
					   G4Colour::Blue(),
					   isCalorimeterSolid,
					   forceAuxEdgeVisible,
					   placePV,
					   doSurfaceCheck
					   );
    if ( verbosityLevel > 0) {
      double zhl         = static_cast<G4Box*>(calorimeterInfo.solid)->GetZHalfLength();
      CLHEP::Hep3Vector const & CalorimeterOffsetInMu2e = calorimeterInfo.centerInMu2e();
      double CalorimeterOffsetInMu2eZ = CalorimeterOffsetInMu2e[CLHEP::Hep3Vector::Z];
      cout << __func__ << " Calorimeter mother center in Mu2e   : " << CalorimeterOffsetInMu2e << endl;
      cout << __func__ << " Calorimeter mother Z extent in Mu2e    : " <<
	CalorimeterOffsetInMu2eZ - zhl << ", " << CalorimeterOffsetInMu2eZ + zhl << endl;
    }

    G4int nRO                   = cal.nROPerCrystal();
    //G4int ncrysR                = cal.nCrystalR();
    //G4int ncrysZ                = cal.nCrystalZ();

    G4double crystalSize        = cal.crystalHalfSize();
    G4double crystalLength      = cal.crystalHalfLength();
    G4double wrapSize           = crystalSize   + cal.wrapperThickness();
    G4double wrapLength         = crystalLength + cal.wrapperThickness() + cal.roHalfThickness();
    G4double shellSize          = wrapSize      + cal.shellThickness();
    G4double shellLength        = wrapLength;


    //-- Create solids for one crystal
    G4Box *crystalShell = new G4Box("CrystalShell",shellLength,shellSize,shellSize);
    G4Box *crystalWrap  = new G4Box("CrystalWrap",wrapLength,wrapSize,wrapSize);
    G4Box *crystal      = new G4Box("Crystal",crystalLength,crystalSize,crystalSize);
    G4Box *crystalRO    = new G4Box("CrystalRO",cal.roHalfThickness(),cal.roHalfSize(),cal.roHalfSize() );


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
    
    if(ccSD) CrystalLog->SetSensitiveDetector(ccSD);
    if(crSD) ROLog->SetSensitiveDetector(crSD);





    //-- Build vanes and absorbers
    const int nvane = cal.nVane();
    const double shieldHalfThickness =  cal.shieldHalfThickness();
    const double neutronAbsorberHalfThickness = cal.neutronAbsorberHalfThickness();
    const double absorberHalfThickness =  shieldHalfThickness + neutronAbsorberHalfThickness;


    VolumeInfo vaneInfo[nvane];
    VolumeInfo absorberInfo[nvane];
    VolumeInfo shieldInfo[nvane];
    VolumeInfo neutronAbsorberInfo[nvane];

    for( int iv=0; iv<nvane; ++iv ) {

      ostringstream nameVane;              nameVane             << "CalorimeterVane_"               << iv;
      ostringstream nameAbsorber;          nameAbsorber         << "CalorimeterAbsorber_"           << iv;
      ostringstream nameShield;            nameShield        << "CalorimeterShield_"          << iv;
      ostringstream nameNeutronAbsorber;   nameNeutronAbsorber  << "CalorimeterNeutronAbsorber_"    << iv;

      const CLHEP::Hep3Vector & size = cal.vane(iv).size();
      double dimVane[3]                = {size.x(), size.y(), size.z() };
      double dimAbsorber[3]            = {size.x(), size.y(), absorberHalfThickness };
      double dimShield[3]              = {size.x(), size.y(), shieldHalfThickness };
      double dimNeutronAbsorber[3]     = {size.x(), size.y(), neutronAbsorberHalfThickness };

      double shift = 0.03;

      G4ThreeVector pVane              = cal.vane(iv).originLocal();
      G4ThreeVector posVane            = G4ThreeVector(pVane.x()-pos.x(), pVane.y()-pos.y(), 0.0);
      G4ThreeVector posAbsorber        = G4ThreeVector(posVane.x(), posVane.y(), posVane.z()-(size.z()+absorberHalfThickness+shift) );
      G4ThreeVector posShield          = G4ThreeVector(0., 0., absorberHalfThickness - shieldHalfThickness);
      G4ThreeVector posNeutronAbsorber = G4ThreeVector(0., 0., absorberHalfThickness - 2*shieldHalfThickness - neutronAbsorberHalfThickness);

      vaneInfo[iv]     = nestBox(nameVane.str(),
				 dimVane,
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

      cout << "Calorimeter Vane position: ("
	   << posVane.x() << "," << posVane.y() << "," << posVane.z()
	   << ")" << endl;

     
      absorberInfo[iv] = nestBox(nameAbsorber.str(),
				 dimAbsorber,
				 fillMaterial,
				 &cal.vane(iv).rotation(),
				 posAbsorber,
				 calorimeterInfo ,
				 iv,
				 isAbsorberBoxVisible,
				 G4Colour::Red(),
				 isAbsorberBoxSolid,
				 forceAuxEdgeVisible,
				 placePV,
				 doSurfaceCheck );

      cout << "Calorimeter Absorber position: ("
	   << posAbsorber.x() << "," << posAbsorber.y() << "," << posAbsorber.z()
	   << ")" << endl;
    
      
      if( shieldHalfThickness > 0.){
	shieldInfo[iv]    = nestBox(nameShield.str(),
				    dimShield,
				    shieldMaterial,
				    0,
				    posShield,
				    absorberInfo[iv] ,
				    iv,
				    isShieldVisible,
				    G4Colour::Green(),
				    isShieldSolid,
				    forceAuxEdgeVisible,
				    placePV,
				    doSurfaceCheck );
      }

      if( neutronAbsorberHalfThickness > 0.){
	neutronAbsorberInfo[iv]    = nestBox(nameNeutronAbsorber.str(),
					     dimNeutronAbsorber,
					     neutronAbsorberMaterial,
					     0,
					     posNeutronAbsorber,
					     absorberInfo[iv] ,
					     iv,
					     isNeutronAbsorberVisible,
					     G4Colour::Magenta(),
					     isNeutronAbsorberSolid,
					     forceAuxEdgeVisible,
					     placePV,
					     doSurfaceCheck );
      }
      
      if ( verbosityLevel > 0 && iv==0) {
	double zhl         = static_cast<G4Box*>(vaneInfo[0].solid)->GetZHalfLength();
	CLHEP::Hep3Vector const & CalorimeterVaneOffsetInMu2e = vaneInfo[0].centerInMu2e();
	double CalorimeterVaneOffsetInMu2eZ = CalorimeterVaneOffsetInMu2e[CLHEP::Hep3Vector::Z];
	cout << __func__ << " Calorimeter mother center in Mu2e   : " << CalorimeterVaneOffsetInMu2e << endl;
	cout << __func__ << " CalorimeterVane   center in Mu2e    : " << CalorimeterVaneOffsetInMu2e << endl;
	cout << __func__ << " CalorimeterVane Z extent in Mu2e    : " <<
	  CalorimeterVaneOffsetInMu2eZ - zhl << ", " << CalorimeterVaneOffsetInMu2eZ + zhl << endl;
	cout<< "zHalfLength = "<< zhl <<endl;
      }





      //-- place crystals inside vanes
      G4int ncrys                 = cal.vane(iv).nCrystals();
      for( int ic=0; ic<ncrys; ++ic ) {

	// IDs
	G4int id       = iv*ncrys + ic;       // Crystal ID
	G4int roidBase = cal.ROBaseByCrystal(id);


	// Have to define a shell / wrapper logical volume for each crystal 
	// to get correct index in CrystalCaloSD
	G4LogicalVolume *thisShellLog(0);
	if (cal.shellThickness() > 0.001){
	  thisShellLog = new G4LogicalVolume(crystalShell, fillMaterial, "ShellLog");
	  thisShellLog->SetVisAttributes(G4VisAttributes::Invisible);
	}

	G4LogicalVolume *thisWrapLog = new G4LogicalVolume(crystalWrap, wrapMaterial, "WrapLog");
	thisWrapLog->SetVisAttributes(G4VisAttributes::Invisible);



	// Position - first run along Z, then along Y, both times in positive direction
	//this is the position of the wrapper, so the x must be zero, not crystalPosition.x()
	CLHEP::Hep3Vector crystalPosition = cal.vane(iv).crystal(ic).position();
	double x = 0;
	double y = crystalPosition.y();
	double z = crystalPosition.z();
	     

	// place a shell only if it has non-zero thickness, or place the wrapper directly
	if (cal.shellThickness()  > 0.001) {
	  new G4PVPlacement(0,G4ThreeVector(x,y,z),thisShellLog,"CrysShellPV",vaneInfo[iv].logical,0,id,doSurfaceCheck);   
	  new G4PVPlacement(0,G4ThreeVector(0.0,0.0,0.0),thisWrapLog,"CrysWrapPV",thisShellLog,0,id,doSurfaceCheck);
	} else {
	  new G4PVPlacement(0,G4ThreeVector(x,y,z),thisWrapLog,"CrysWrapPV",vaneInfo[iv].logical,0,id,doSurfaceCheck);   	      
	}

	// -- place crystal inside warp
	new G4PVPlacement(0,G4ThreeVector(cal.roHalfThickness(),0.0,0.0),CrystalLog,"CrysPV",thisWrapLog,0,id,doSurfaceCheck);

	// -- add the readout
	if (nRO==1) 
	  new G4PVPlacement(0,G4ThreeVector(-cal.crystalHalfLength(),0,0),ROLog,"CrysROPV_0",thisWrapLog,0,roidBase,doSurfaceCheck);

	if (nRO==2) { 
	  new G4PVPlacement(0,G4ThreeVector(-cal.crystalHalfLength(),-0.5*cal.crystalHalfSize(),0.0),ROLog,"CrysROPV_0",thisWrapLog,0,roidBase,doSurfaceCheck);
	  new G4PVPlacement(0,G4ThreeVector(-cal.crystalHalfLength(), 0.5*cal.crystalHalfSize(),0.0),ROLog,"CrysROPV_1",thisWrapLog,0,roidBase+1,doSurfaceCheck);
	}

	if (nRO==4) { 
	  G4double cHS = -0.5*cal.crystalHalfSize();
	  new G4PVPlacement(0,G4ThreeVector(-cal.crystalHalfLength(),-cHS,-cHS),ROLog,"CrysROPV_0",thisWrapLog,0,roidBase,doSurfaceCheck);
	  new G4PVPlacement(0,G4ThreeVector(-cal.crystalHalfLength(),-cHS,cHS), ROLog,"CrysROPV_1",thisWrapLog,0,roidBase+1,doSurfaceCheck);
	  new G4PVPlacement(0,G4ThreeVector(-cal.crystalHalfLength(), cHS,-cHS),ROLog,"CrysR0PV_2",thisWrapLog,0,roidBase+2,doSurfaceCheck);
	  new G4PVPlacement(0,G4ThreeVector(-cal.crystalHalfLength(), cHS,cHS), ROLog,"CrysROPV_3",thisWrapLog,0,roidBase+3,doSurfaceCheck);
	}


      }//end crystal loop



    }//end loop over vanes


    return calorimeterInfo;

  }


} // end namespace mu2e
