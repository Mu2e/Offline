//
// Free function to create the calorimeter.
//
// $Id: constructCalorimeter.cc,v 1.14 2011/05/17 15:36:01 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:01 $
//
// Original author Ivan Logashenko
// 
// Notes
// 1) The argument zOff is the zlocation of the center of the mother volume,
//    as mesaured in the mu2e coordinate system.
// 2) The vanes are placed directly into DS3.  We did not make a mother volume for them.

#include <iostream>

// Mu2e includes.
#include "Mu2eG4/inc/constructCalorimeter.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
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

  void constructCalorimeter( VolumeInfo const &  mother,
                             double              zOffset,
                             SimpleConfig const& config ){

    int const verbosityLevel = config.get<int>("calorimeter.verbosityLevel",0);

    // A helper class for parsing the config file.
    MaterialFinder materialFinder(config);

    Calorimeter const & cal = *(GeomHandle<Calorimeter>());

    // Read parameters from config file
    bool isVaneBoxVisible     = config.get<bool>("calorimeter.vaneBoxVisible",true);
    bool isVaneBoxSolid       = config.get<bool>("calorimeter.vaneBoxSolid",true);
    bool isCrystalVisible     = config.get<bool>("calorimeter.crystalVisible",false);
    bool isCrystalSolid       = config.get<bool>("calorimeter.crystalSolid",true);
    bool forceAuxEdgeVisible  = config.get<bool>("g4.forceAuxEdgeVisible",false);
    G4bool doSurfaceCheck     = config.get<bool>("g4.doSurfaceCheck",false);
    bool const placePV = true;
    int roSide = config.get<int>("calorimeter.crystalReadoutSide",1);
    if( roSide>=0 ) roSide=1; else roSide=-1;

    G4Material* fillMaterial = materialFinder.get("calorimeter.calorimeterFillMaterial");

    const int nvane = cal.nVane();
    VolumeInfo vaneInfo[nvane];

    G4ThreeVector pcalo = cal.getOrigin();

    for( int i=0; i<nvane; ++i ) {

      G4ThreeVector pvane = cal.getVane(i).getOriginLocal();
      G4ThreeVector pos  = G4ThreeVector(pvane.x(), pvane.y(), pcalo.z()+zOffset);

      const CLHEP::Hep3Vector & size = cal.getVane(i).getSize();

      double dim[3] = { size.x(), size.y(), size.z() };

      cout << "Calorimeter Vane position: (" 
           << pos.x() << "," << pos.y() << "," << pos.z() 
           << ")" << endl;

      ostringstream name;
      name << "CalorimeterVane_" << i;

      vaneInfo[i] = nestBox(name.str(),
                            dim,
                            fillMaterial,
                            cal.getVane(i).getRotation(),
                            pos,
                            mother,
                            i,
                            isVaneBoxVisible,
                            G4Colour::Yellow(),
                            isVaneBoxSolid,
                            forceAuxEdgeVisible,
                            placePV,
                            doSurfaceCheck );

    }

    if ( verbosityLevel > 0) {
      double zhl         = static_cast<G4Box*>(vaneInfo[0].solid)->GetZHalfLength();
      CLHEP::Hep3Vector const & CalorimeterVaneOffsetInMu2e = vaneInfo[0].centerInMu2e();
      double CalorimeterVaneOffsetInMu2eZ = CalorimeterVaneOffsetInMu2e[CLHEP::Hep3Vector::Z];
      cout << __func__ << " Calorimeter mother center in Mu2e   : " << CalorimeterVaneOffsetInMu2e << endl;
      cout << __func__ << " CalorimeterVane   center in Mu2e    : " << CalorimeterVaneOffsetInMu2e << endl;
      cout << __func__ << " CalorimeterVane Z extent in Mu2e    : " << 
        CalorimeterVaneOffsetInMu2eZ - zhl << ", " << CalorimeterVaneOffsetInMu2eZ + zhl << endl;
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

    // Create logical volumes - we are going to reuse these
        
    G4LogicalVolume * l_wrap  = new G4LogicalVolume( wrap,  wrapMaterial, "l_CrystalWrap"); 
    l_wrap->SetVisAttributes(G4VisAttributes::Invisible);

    G4LogicalVolume * l_crys  = new G4LogicalVolume( crys,  crysMaterial, "l_Crystal"); 
    l_crys->SetVisAttributes(new G4VisAttributes(false,G4Colour::Green()));

    G4VSensitiveDetector* ccSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::CaloCrystal());
    l_crys->SetSensitiveDetector(ccSD);

    G4LogicalVolume * l_ro = new G4LogicalVolume( ro,    readMaterial,  "l_CrystalRO" ); 
    if(!isCrystalVisible) {
      l_ro->SetVisAttributes(G4VisAttributes::Invisible);
    } else {
      G4VisAttributes* ro_visAtt = new G4VisAttributes(isCrystalVisible, G4Color::Red());
      ro_visAtt->SetForceSolid(isCrystalSolid);
      ro_visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
      l_ro->SetVisAttributes(ro_visAtt);
    }

    G4VSensitiveDetector* crSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::CaloReadout());
    l_ro->SetSensitiveDetector(crSD);


    //
    // Create crystals in the loop, one at a time. I do it this way, not using
    // replica or reusing the same logical volume to avoid problem with volume
    // index in physHelper. 
    //

    double step = cal.crystalHalfSize()*2.0;

    for( int iv=0; iv<nvane; ++iv ) {
      for( int ic=0; ic<ncrys; ++ic ) {

        // IDs
        int id     = iv*ncrys + ic;       // Crystal ID
        int roid   = nro*(iv*ncrys + ic); // Readout ID

        ostringstream lname; lname << "l_CrystalShell" << id;

        G4LogicalVolume * l_shell = new G4LogicalVolume( shell, fillMaterial, lname.str()); 
        if(!isCrystalVisible) {
          l_shell->SetVisAttributes(G4VisAttributes::Invisible);
        } else {
          G4VisAttributes* shell_visAtt = new G4VisAttributes(isCrystalVisible, G4Color::Blue());
          shell_visAtt->SetForceSolid(isCrystalSolid);
          shell_visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
          l_shell->SetVisAttributes(shell_visAtt);
        }

        ostringstream wname; wname << "l_CrystalWrap" << id;

        G4LogicalVolume * l_wrap  = new G4LogicalVolume( wrap,  wrapMaterial, wname.str()); 
        l_wrap->SetVisAttributes(G4VisAttributes::Invisible);

        // -- place crystal inside wrap
        new G4PVPlacement(0,G4ThreeVector(0.0,0.0,0.0),l_crys,
                          "p_Crystal",l_wrap,0,id,doSurfaceCheck);

        // -- place wrap inside shell
        // p_CrystalShell not p_CrystalWrap ???
        new G4PVPlacement(0,G4ThreeVector(roSide*cal.roHalfThickness(),0.0,0.0),l_wrap,
                          "p_CrystalShell",l_shell,0,id,doSurfaceCheck);

        // -- add readouts
        for( int i=0; i<nro; ++i ) {
          if( nro==1 ) {
            new G4PVPlacement(0,G4ThreeVector(-roSide*cal.crystalHalfLength(),0.0,0.0),
                              l_ro,"p_CrystalRO",l_shell,0,roid+i,doSurfaceCheck);
          } else if( nro==2 ) {
            new G4PVPlacement(0,G4ThreeVector(-roSide*cal.crystalHalfLength(),(i-0.5)*cal.crystalHalfSize(),0.0),
                              l_ro,"p_CrystalRO",l_shell,0,roid+i,doSurfaceCheck);
          } else if( nro==4 ) {
            new G4PVPlacement(0,G4ThreeVector(-roSide*cal.crystalHalfLength(),
                                              (i/2-0.5)*cal.crystalHalfSize(),
                                              (i%2-0.5)*cal.crystalHalfSize()),
                              l_ro,"p_CrystalRO",l_shell,0,roid+i,doSurfaceCheck);
          }
        }

        // Place crystal shell for each crystal in the vane. 
        
        // Position - first run along Z, then along Y, both times in positive direction
        double x = 0.0;
        double y = 0.5*step*(2*(ic/ncrysZ)-ncrysR+1);
        double z = 0.5*step*(2*(ic%ncrysZ)-ncrysZ+1);

        // Create volumes 
        
        new G4PVPlacement(0,G4ThreeVector(x,y,z),l_shell,"Crystal"
                          ,vaneInfo[iv].logical,0,id,doSurfaceCheck);
      }
    }
    
  }
  
} // end namespace mu2e
