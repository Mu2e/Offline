#ifndef Mu2eWorld_H
#define Mu2eWorld_H 1
//
// Construct the Mu2e G4 world and serve information about that world.
//
// $Id: Mu2eWorld.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include <string>
#include <memory>

// Forward references.
class G4Material;
class G4UniformMagField;
class G4Mag_UsualEqRhs;
class G4ExactHelixStepper;
class G4ChordFinder;
class G4FieldManager;
class G4UserLimits;

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Colour.hh"

#include "Mu2eG4/inc/WorldInfo.hh"
#include "Mu2eG4/inc/VolumeInfo.hh"

#include "G4ThreeVector.hh"

namespace mu2e {

  // Forward references within mu2e namespace.
  class SimpleConfig;

  class Mu2eWorld {
  public:
    
    Mu2eWorld();
    ~Mu2eWorld();

    // Construct everything.
    WorldInfo const* construct();

    WorldInfo const& getWorldInfo() const { return _info; }

    G4ThreeVector const& getCosmicReferencePoint() const{
      return _cosmicReferencePoint;
    }

    G4ThreeVector const& getMu2eOrigin() const{
      return _mu2eOrigin;
    }

    G4ThreeVector const& getMu2eDetectorOrigin() const{
      return _mu2eDetectorOrigin;
    }

  private:
    
    void constructWorld( SimpleConfig const& );
    VolumeInfo constructLTracker( G4LogicalVolume* mother, double zOff );
    void constructTestWorld();

    // The world coordinates of the center of the cosmic ray reference plane.
    G4ThreeVector _cosmicReferencePoint;

    // The world coordinates of the origin of the Mu2e coordinate system.
    G4ThreeVector _mu2eOrigin;

    // World coorindates of the reference point for building the detectors:
    //  - on axis on the DS with z=12,000. in the Mu2e coordinate system.
    G4ThreeVector _mu2eDetectorOrigin;

    // Information about the world that can be passed to others.
    WorldInfo _info;

    // Utility functions.
    void Mu2eWorld::setUnits( std::vector<double>& V, G4double unit );
    G4Material* Mu2eWorld::getMaterial( std::string const& name );

    SimpleConfig const* _config;

    //Information for a few upper level volumes.
    //    G4Box* _worldSolid;
    
    // Logical volumes
    //G4LogicalVolume* _worldLog;
    
    // Physical volumes
    //G4VPhysicalVolume* _worldPhys;

    std::auto_ptr<G4UniformMagField>   _detSolBField;
    std::auto_ptr<G4Mag_UsualEqRhs>    _usualRHS;
    std::auto_ptr<G4ExactHelixStepper> _exactHelix;
    std::auto_ptr<G4ChordFinder>       _chordFinder;
    std::auto_ptr<G4FieldManager>      _fieldMgr;
    std::auto_ptr<G4UserLimits>        _stepLimit;


  };

} // end namespace mu2e
#endif

