#ifndef Mu2eWorld_H
#define Mu2eWorld_H 1
//
// Construct the Mu2e G4 world and serve information about that world.
//
// $Id: Mu2eWorld.hh,v 1.12 2010/05/18 22:34:46 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/05/18 22:34:46 $
//
// Original author Rob Kutschke
//

#include <string>
#include <memory>

// Forward references.
class G4Material;
class DSField;
class G4UniformMagField;
class G4Mag_UsualEqRhs;
class G4ExactHelixStepper;
class G4ChordFinder;
class G4FieldManager;
class G4UserLimits;
class G4HelixSimpleRunge;
class G4ClassicalRK4;
class G4CashKarpRKF45;
class G4ImplicitEuler;
class G4ExplicitEuler;
class G4MagneticField;

//
//G4 includes 
#include "G4String.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

//
// Mu2e includes
#include "Mu2eG4/inc/WorldInfo.hh"
#include "Mu2eG4/inc/VolumeInfo.hh"



class G4AssemblyVolume;

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

    G4ThreeVector const& getPrimaryProtonGunOrigin() const{
      return _primaryProtonGunOrigin;
    }
 
    G4RotationMatrix const& getPrimaryProtonGunRotation() const{
      return _primaryProtonGunRotation;
    }


  private:
    
    void constructWorld( SimpleConfig const& );

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

    SimpleConfig const* _config;

    // Location in G4 world coordinates of the reference point for the Primary Proton Gun
    G4ThreeVector _primaryProtonGunOrigin;
    G4RotationMatrix _primaryProtonGunRotation;

    //keep these for future use
    std::auto_ptr<G4UniformMagField>   _detSolUpstreamBField;
    std::auto_ptr<G4UniformMagField>   _detSolDownstreamBField;

    //enum to pick field form
    enum detSolFieldChoice {detSolFullField,detSolUpVaryingDownConstant,detSolUpConstantDownConstant};

    //
    //if you want to use constant field
    std::auto_ptr<G4UniformMagField>   _detSolUpstreamConstantBField;
    std::auto_ptr<G4UniformMagField>   _detSolDownstreamConstantBField;
    //
    //varying field
    std::auto_ptr<G4MagneticField>   _detSolUpstreamVaryingBField;
    std::auto_ptr<G4MagneticField>   _detSolDownstreamVaryingBField;


    //need these in both cases
    std::auto_ptr<G4Mag_UsualEqRhs>    _usualUpstreamRHS;
    std::auto_ptr<G4ExactHelixStepper> _exactUpstreamHelix;
    std::auto_ptr<G4ClassicalRK4 >     _rungeUpstreamHelix;
    std::auto_ptr<G4CashKarpRKF45 >     _rungeCKUpstreamHelix;
    std::auto_ptr<G4ImplicitEuler >     _rungeIEUpstreamHelix;
    std::auto_ptr<G4ExplicitEuler >     _rungeEEUpstreamHelix;
    std::auto_ptr<G4ChordFinder>       _chordUpstreamFinder;
    std::auto_ptr<G4FieldManager>      _fieldUpstreamMgr;
    std::auto_ptr<G4Mag_UsualEqRhs>    _usualDownstreamRHS;
    std::auto_ptr<G4ExactHelixStepper> _exactDownstreamHelix;
    std::auto_ptr<G4ClassicalRK4>     _rungeDownstreamHelix;
    std::auto_ptr<G4ExplicitEuler>     _rungeEEDownstreamHelix;
    std::auto_ptr<G4ChordFinder>       _chordDownstreamFinder;
    std::auto_ptr<G4FieldManager>      _fieldDownstreamMgr;
    std::auto_ptr<G4UserLimits>        _stepLimit;
    std::auto_ptr<G4UserLimits>        _stepUpstreamLimit;
    std::auto_ptr<G4UserLimits>        _stepDownstreamLimit;


  };

} // end namespace mu2e
#endif
