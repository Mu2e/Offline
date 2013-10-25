//
// Construct the Mu2e detector with the Mu2e G4 world.
//
// $Id: DetectorConstruction.cc,v 1.6 2013/10/25 18:37:07 genser Exp $
// $Author: genser $
// $Date: 2013/10/25 18:37:07 $
//
// Original author Rob Kutschke
//

#include <iostream>
#include <sstream>

#include "Mu2eG4/inc/DetectorConstruction.hh"

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// Geant4 includes
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ChordFinder.hh"

#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

using namespace std;

namespace mu2e {

  DetectorConstruction::DetectorConstruction()
    :_magField(0),
     _stepLimit(0),
     _usualRHS(0),
     _exactHelix(0),
     _chordFinder(0){
  }


  DetectorConstruction::~DetectorConstruction(){

    // For most of the objects that are new'ed in the construct method,
    // G4 takes ownship of the object and will delete it when finished.
    //
    // The exceptions are:

    delete _magField;
    delete _stepLimit;
    delete _usualRHS;
    delete _exactHelix;
    delete _chordFinder;

  }


  G4VPhysicalVolume* DetectorConstruction::Construct(){

    // Some made up parameters to describe vacuum.
    G4double density     = 1.e-10*CLHEP::g/CLHEP::cm3;
    G4double pressure    = 3.e-18*CLHEP::pascal;
    G4double temperature = 2.73*CLHEP::kelvin;
    G4Material* vacuum =
      new G4Material( "Vacuum", 1., 1.01*CLHEP::g/CLHEP::mole,
                      density, kStateGas, temperature, pressure);

    // Print all the materials defined.
    G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;

    // Half dimensions of the world.  Roughly the size of the Detector Solenoid.
    _xwHalf = 1000.*CLHEP::mm;
    _ywHalf = 1000.*CLHEP::mm;
    _zwHalf = 3000.*CLHEP::mm;

    // Thickness in z of the reference boxes.
    const G4double zHalfThick(0.010*CLHEP::mm);

    // Number of reference volumes to make;
    const G4int nRef(11);

    // Offset of the first and last reference volumes from the edge of the world
    const double zoff(_zwHalf*0.05);

    // Color for visualizing the reference boxes
    const G4VisAttributes visRefBox(true,G4Color::Red());

    // Magnetic Vector.
    const G4ThreeVector bfield(0.,0.,1.*CLHEP::tesla);

    // Compute some derived parameters.

    // Z position of the negative z edge of the first reference box.
    const G4double z0 = -_zwHalf+zoff;

    // Starting position for tracks.
    _zStart = z0-zHalfThick;

    // Distance in z between reference boxes.
    // N planes = (N-1) spaces between planes.
    _dz = 2.*(_zwHalf-zoff)/(nRef-1);

    // End of the parameter definitions.  Start making things.

    // Make the global magnetic field.
    G4UniformMagField* _magField = new G4UniformMagField(bfield);

    // Tell the global field manager to use this field.
    G4FieldManager* fieldMgr= G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(_magField);

    // Define propagation code that does an exact helix in a uniform magnetic field.
    G4double stepMinimum(1.0e-2*CLHEP::mm);
    _usualRHS    = new G4Mag_UsualEqRhs( _magField );
    _exactHelix  = new G4ExactHelixStepper(_usualRHS);
    _chordFinder = new G4ChordFinder( _magField, stepMinimum, _exactHelix );

    // Tell the field manager to use the exact propagation code.
    fieldMgr->SetChordFinder(_chordFinder);

    // Make the world.

    G4Box* solidWorld = new G4Box( "world", _xwHalf, _ywHalf, _zwHalf );
    G4LogicalVolume* logicWorld = new G4LogicalVolume( solidWorld, vacuum, "World", 0, 0, 0 );
    G4VPhysicalVolume* physiWorld =
      new G4PVPlacement(0,               // no rotation
                        G4ThreeVector(), // at (0,0,0)
                        logicWorld,      // its logical volume
                        "World",         // its name
                        0,               // its mother  volume
                        false,           // no boolean operations
                        0);              // copy number

    // Divide the work into two hemispheres.




    // Give each box a name.
    G4String sname = "Ref_Solid";
    G4String lname = "Ref_Log";
    G4String pname = "Ref_Phys";


    // Make a solid.
    G4Box* solid = new G4Box( sname, _xwHalf, _ywHalf, zHalfThick);

    // Make a logical volume and set its color.
    G4LogicalVolume* logical  = new G4LogicalVolume( solid,
                                                     vacuum,
                                                     lname);
    logical->SetVisAttributes(visRefBox);



    // Create and position the reference boxes.
    for ( int i=0; i<nRef; ++i ){

      // Position of this reference box.
      double z = z0 +i*_dz;
      G4ThreeVector position( 0., 0., z);

      // Place the volume.
      new G4PVPlacement( 0,
                         position,
                         logical,
                         pname,
                         logicWorld,
                         false,
                         i
                         );
    }

    // Set a maximum set length.  This does not seem to work.
    _stepLimit = new G4UserLimits(5.*CLHEP::mm);
    logicWorld->SetUserLimits(_stepLimit);

    return physiWorld;

  }

} // end namespace mu2e
