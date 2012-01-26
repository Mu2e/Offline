#ifndef Mu2eG4_Mu2eWorld_hh
#define Mu2eG4_Mu2eWorld_hh
//
// Construct the Mu2e G4 world and serve information about that world.
//
// $Id: Mu2eWorld.hh,v 1.39 2012/01/26 21:57:12 genser Exp $
// $Author: genser $
// $Date: 2012/01/26 21:57:12 $
//
// Original author Rob Kutschke
//
// Notes
// 1) The variable _volumeInfoList, holds some pointers to information
//    about volumes.  It also holds some position information.  The
//    data member centerInWorld is not always meaningful.  If you follow
//    the trail from a given volume back to the world volume, if you
//    do not find any rotations, then this data member is meaningful.
//    If you do find a rotation, then the data member is not meaningful.
//    This data member will only be filled for some of the upper level
//    volumes.   It won't be filled for straws and crystals.  It's purpose
//    is to make it easier to break up one giant method into many smaller
//    ones, by allowing code to look up its mother volume by name.
//    The bottom line is that this is not a fully general facility and
//    must be used with care.

#include <string>
#include <memory>
#include <vector>
#include <map>

// Forward references.
class G4Material;
class G4Mag_UsualEqRhs;
class G4UserLimits;

// Mu2e includes
#include "G4Helper/inc/VolumeInfo.hh"
#include "Mu2eG4/inc/FieldMgr.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "TrackerGeom/inc/TubsParams.hh"

//G4 includes
#include "G4String.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"

namespace mu2e {

  // Forward references within mu2e namespace.
  class SimpleConfig;

  class Mu2eWorld {
  public:

    Mu2eWorld();
    ~Mu2eWorld();

    // Construct everything.
    // The non-const return type is eventually required 
    // by G4VUserDetectorConstruction::Construct();
    G4VPhysicalVolume * construct();

  private:

    // Do all of the work.
    G4VPhysicalVolume * constructWorld();

    // Break the big task into many smaller ones.
    VolumeInfo constructTracker();
    VolumeInfo constructTarget();
    void constructCal();
    void constructMagnetYoke();
    void constructBFieldAndManagers();
    void constructStepLimiters();
    void constructITStepLimiters();

    void instantiateSensitiveDetectors();

    // Utility functions.
    void setUnits( std::vector<double>& V, G4double unit );

    // Stash a pointer to the config object so that all methods can get at it easily.
    SimpleConfig const* _config;

    // Models of the DS magnetic field:
    // 0 - whole DS uses the field map.
    // 1 - upstream uses the full field map; downstream uses a uniform field.
    // 2 - whole DS uses a uniform field.
    enum DSFieldModel { dsModelFull, dsModelSplit, dsModelUniform};

    // Field managers for the different regions of magnetic field.
    // These have a lifetime equal to that of the G4 geometry.
    std::auto_ptr<FieldMgr> _dsUniform;
    std::auto_ptr<FieldMgr> _dsGradient;

    // Access to the G4HelperService.
    G4Helper * _helper;

    int  _verbosityLevel;

  };

} // end namespace mu2e
#endif /* Mu2eG4_Mu2eWorld_hh */
