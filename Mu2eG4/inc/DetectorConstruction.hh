#ifndef Mu2eG4_DetectorConstruction_hh
#define Mu2eG4_DetectorConstruction_hh
//
// Construct the Mu2e detector with the Mu2e G4 world.
//
// $Id: DetectorConstruction.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4UniformMagField.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4UserLimits;
class G4Mag_UsualEqRhs;
class G4ExactHelixStepper;
class G4ChordFinder;

namespace mu2e {

  class DetectorConstruction : public G4VUserDetectorConstruction{
  public:
    
    DetectorConstruction();
    ~DetectorConstruction();
    
  public:
    
    G4VPhysicalVolume* Construct();

    G4double getZStart() const { return _zStart;}
    G4double getDz() const { return _dz;}
    G4double getXHalf() const { return _xwHalf;}
    G4double getYHalf() const { return _ywHalf;}
    G4double getZHalf() const { return _zwHalf;}

  private:

    // Pointer to the magnetic field 
    G4UniformMagField* _magField;
    G4UserLimits*      _stepLimit;

    // Stuff needed to use the exact field code.
    G4Mag_UsualEqRhs    *_usualRHS;
    G4ExactHelixStepper *_exactHelix;
    G4ChordFinder       *_chordFinder;
    
    // Half lengths of the world.
    G4double _xwHalf;
    G4double _ywHalf;
    G4double _zwHalf;

    // Starting z for tracks.
    G4double _zStart;

    // Delta z between reference planes.
    G4double _dz;

  };

}  // end namespace mu2e
#endif /* Mu2eG4_DetectorConstruction_hh */
