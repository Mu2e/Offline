//
// Free function to create the calorimeter.
//
// $Id: constructCalorimeter.cc,v 1.8 2010/12/02 17:49:21 genser Exp $
// $Author: genser $
// $Date: 2010/12/02 17:49:21 $
//
// Original author Rob Kutschke
// 
// Notes
// 1) The argument zOff is the zlocation of the center of the mother volume,
//    as mesaured in the mu2e coordinate system.

#include <iostream>

// Mu2e includes.
#include "Mu2eG4/inc/constructCalorimeter.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/Vane.hh"
#include "Mu2eG4/inc/CaloCrystalSD.hh"
#include "Mu2eG4/inc/CaloReadoutSD.hh"

// G4 includes
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"

using namespace std;

namespace mu2e {

  void constructCalorimeter( G4LogicalVolume*    mother,
                                   double              zOffset,
                                   SimpleConfig const& config ){
    // A helper class for parsing the config file.
    MaterialFinder materialFinder(config);

    Calorimeter const & cal = *(GeomHandle<Calorimeter>());

    // Read parameters from config file
    bool isVisible = config.getBool("calorimeter.visible",true);
    bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    bool isSolid   = config.getBool("calorimeter.solid",true);
    G4bool doSurfaceCheck = config.getBool("g4.doSurfaceCheck",false);
    bool const placePV = true;

    G4Material* fillMaterial = materialFinder.get("calorimeter.calorimeterFillMaterial");
    //G4Material* fillMaterial = materialFinder.get("calorimeter.crystalMaterial");

    // Create vanes. Do not create mother volume for calorimeter - 
    // add vanes directly to DS3 (mother).

    const int nvane = cal.nVane();

    VolumeInfo vaneInfo[nvane];

    G4ThreeVector pcalo = cal.getOrigin();

    for( int i=0; i<nvane; ++i ) {

      G4ThreeVector pvane = cal.getVane(i).getOriginLocal();
      G4ThreeVector pos  = G4ThreeVector(pvane.x(), pvane.y(), pcalo.z()+zOffset);

      const CLHEP::Hep3Vector & size = cal.getVane(i).getSize();

      double dim[3] = { size.x(), size.y(), size.z() };

      cout << "Calorimeter Vane position: (" << pos.x() << "," << pos.y() << "," << pos.z() << ")" << endl;

      ostringstream name;
      name << "CalorimeterVane_" << i;

      vaneInfo[i] = nestBox(name.str(), 
                            dim,
                            fillMaterial,
			    cal.getVane(i).getRotation(), 
                            pos,
			    mother, 
                            i,
                            isVisible,
			    G4Colour::Yellow(), 
                            isSolid,
                            forceAuxEdgeVisible,
                            placePV,
                            doSurfaceCheck );

    }

    // Create materials for crystals 

    G4Material* crysMaterial = materialFinder.get("calorimeter.crystalMaterial");
    G4Material* wrapMaterial = materialFinder.get("calorimeter.crystalWrapper");
    G4Material* readMaterial = materialFinder.get("calorimeter.crystalReadoutMaterial");

    // Create solids for one crystal

    G4Box * shell = new G4Box("CrystalShell", 
			      cal.crystalHalfLength()+cal.roHalfThickness(), 
			      cal.crystalHalfSize(),
			      cal.crystalHalfSize() );
    G4Box * wrap  = new G4Box("CrystalWrap", 
			      cal.crystalHalfLength(), 
			      cal.crystalHalfSize(),
			      cal.crystalHalfSize() );
    G4Box * crys  = new G4Box("Crystal", 
			      cal.crystalHalfLength()-2.0*cal.wrapperHalfThickness(), 
			      cal.crystalHalfSize()-2.0*cal.wrapperHalfThickness(), 
			      cal.crystalHalfSize()-2.0*cal.wrapperHalfThickness() );
    G4Box * ro    = new G4Box("CrystalRO", 
			      cal.roHalfThickness(), 
			      cal.roHalfSize(),
			      cal.roHalfSize() );
    
    int nro    = cal.nROPerCrystal();
    int ncrys  = cal.nCrystalPerVane();
    int ncrysR = cal.nCrystalR();
    int ncrysZ = cal.nCrystalZ();

    // Create sensitive detectors

    G4SDManager* SDman    = G4SDManager::GetSDMpointer();

    CaloCrystalSD* crysSD = new CaloCrystalSD("CaloCrystal", config);
    SDman->AddNewDetector(crysSD);

    CaloReadoutSD* roSD = new CaloReadoutSD("CaloReadout", config);
    SDman->AddNewDetector(roSD);

    // Create logical volumes
	
    G4LogicalVolume * l_shell = new G4LogicalVolume( shell, fillMaterial, "l_CrystalShell"); 
    if(!isVisible) {
      l_shell->SetVisAttributes(G4VisAttributes::Invisible);
    } else {
      G4VisAttributes* shell_visAtt = new G4VisAttributes(isVisible, G4Color::Blue());
      shell_visAtt->SetForceSolid(isSolid);
      shell_visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
      l_shell->SetVisAttributes(shell_visAtt);
    }

    G4LogicalVolume * l_wrap  = new G4LogicalVolume( wrap,  wrapMaterial, "l_CrystalWrap"); 
    l_wrap->SetVisAttributes(G4VisAttributes::Invisible);

    G4LogicalVolume * l_crys  = new G4LogicalVolume( crys,  crysMaterial, "l_Crystal"); 
    l_crys->SetVisAttributes(new G4VisAttributes(false,G4Colour::Green()));
    l_crys->SetSensitiveDetector(crysSD);

    G4LogicalVolume * l_ro = new G4LogicalVolume( ro,    readMaterial,  "l_CrystalRO" ); 
    if(!isVisible) {
      l_ro->SetVisAttributes(G4VisAttributes::Invisible);
    } else {
      G4VisAttributes* ro_visAtt = new G4VisAttributes(isVisible, G4Color::Red());
      ro_visAtt->SetForceSolid(isSolid);
      ro_visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
      l_ro->SetVisAttributes(ro_visAtt);
    }
    l_ro->SetSensitiveDetector(roSD);

    //
    // Create single crystal
    //
    // -- place crystal inside wrap
    new G4PVPlacement(0,G4ThreeVector(0.0,0.0,0.0),l_crys,
                      "p_Crystal",l_wrap,0,0,doSurfaceCheck);
    // -- place wrap inside shell
    // p_CrystalShell not p_CrystalWrap ???
    new G4PVPlacement(0,G4ThreeVector(-cal.roHalfThickness(),0.0,0.0),l_wrap,
		      "p_CrystalShell",l_shell,0,0,doSurfaceCheck);
    // -- add readouts
    for( int i=0; i<nro; ++i ) {
      ostringstream pname; pname << "p_CrystalRO" << i;
      if( nro==1 ) {
	new G4PVPlacement(0,G4ThreeVector(cal.crystalHalfLength(),0.0,0.0),
			  l_ro,pname.str(),l_shell,0,i,doSurfaceCheck);
      } else if( nro==2 ) {
	new G4PVPlacement(0,G4ThreeVector(cal.crystalHalfLength(),(i-0.5)*cal.crystalHalfSize(),0.0),
			  l_ro,pname.str(),l_shell,0,i,doSurfaceCheck);
      } else if( nro==4 ) {
	new G4PVPlacement(0,G4ThreeVector(cal.crystalHalfLength(),
					  (i/2-0.5)*cal.crystalHalfSize(),
					  (i%2-0.5)*cal.crystalHalfSize()),
			  l_ro,pname.str(),l_shell,0,i,doSurfaceCheck);
      }
    }

    // Place crystal shell for each crystal. If neccessary, this code can be 
    // rewritten to use Replica

    double step = cal.crystalHalfSize()*2.0;

    for( int iv=0; iv<nvane; ++iv ) {
      for( int ic=0; ic<ncrys; ++ic ) {

	// IDs
	int id   = iv*ncrys + ic; // Crystal ID
	
	// Position - first run along Z, then along Y, both times in positive direction
	double x = 0.0;
	double y = 0.5*step*(2*(ic/ncrysZ)-ncrysR+1);
	double z = 0.5*step*(2*(ic%ncrysZ)-ncrysZ+1);

	ostringstream name;
	name << "Crystal" << id;
	
	// Create volumes 

	new G4PVPlacement(0,G4ThreeVector(x,y,z),l_shell,name.str()
                          ,vaneInfo[iv].logical,0,id,doSurfaceCheck);
      }
    }

  }

} // end namespace mu2e
