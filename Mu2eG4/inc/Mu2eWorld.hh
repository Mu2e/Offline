#ifndef Mu2eG4_Mu2eWorld_hh
#define Mu2eG4_Mu2eWorld_hh
//
// Construct the Mu2e G4 world and serve information about that world.
// Note that the class inherits from Mu2eUniverse now
//
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
#include <list>

// Forward references.
class G4Material;
class G4Mag_UsualEqRhs;
class G4UserLimits;

// Mu2e includes
#include "Mu2eG4/inc/Mu2eUniverse.hh"
#include "Mu2eG4/inc/Mu2eG4Config.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "Mu2eG4/inc/FieldMgr.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"


//G4 includes
#include "G4String.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4LogicalVolumeStore;
class G4VSensitiveDetector;

namespace mu2e {

  // Forward references within mu2e namespace.
  class SimpleConfig;
  class SensitiveDetectorHelper;

  class Mu2eWorld : public Mu2eUniverse {
  public:

    Mu2eWorld(const Mu2eG4Config::Top& conf,
              SensitiveDetectorHelper *sdHelper/*no ownership passing*/
              );

    // Construct everything.
    // The non-const return type is eventually required
    // by G4VUserDetectorConstruction::Construct();
    virtual G4VPhysicalVolume * construct() override;

    virtual void constructSDandField() override;

  private:

    typedef std::list<G4LogicalVolume*> LVList;
    typedef std::map<G4VSensitiveDetector*, LVList> DetList;

    // Do all of the work.
    G4VPhysicalVolume * constructWorld();

    // Break the big task into many smaller ones.
    VolumeInfo constructTracker();
    VolumeInfo constructTarget();
    VolumeInfo constructCal();
    void constructMagnetYoke();
    void constructBFieldAndManagers();
    void constructStepLimiters();
    void constructITStepLimiters();

    void instantiateSensitiveDetectors();

    void stepLimiterHelper(const std::string &regexp, G4UserLimits* stepLimit);
    void setStepLimitToAllSuchVolumes(const G4String& vn,
                                      G4UserLimits* const stepLimit,
                                      const G4LogicalVolumeStore* const lvs,
                                      int verbosityLevel);

    // Field managers for the different regions of magnetic field.
    // These have a lifetime equal to that of the G4 geometry.
    std::unique_ptr<FieldMgr> _dsUniform;
    std::unique_ptr<FieldMgr> _dsGradient;

    SensitiveDetectorHelper *sdHelper_; // Non-owning

    Mu2eG4Config::Top conf_;

    // Values of the following variables are taken from either
    // ParameterSet or SimpleConfig, depending on the constructor
    // called.

    // _verbosityLevel in the base class
    bool activeWr_Wl_SD_;
    bool writeGDML_;
    std::string gdmlFileName_;
    std::string g4stepperName_;
    double g4epsilonMin_;
    double g4epsilonMax_;
    double g4DeltaOneStep_;
    double g4DeltaIntersection_;
    double g4DeltaChord_;
    double g4StepMinimum_;
    int    g4MaxIntSteps_;
    double bfieldMaxStep_;
    double strawGasMaxStep_;
    bool limitStepInAllVolumes_;
    bool useEmOption4InTracker_;

    //returned from constructPS
    G4LogicalVolume* psVacuumLogical_;

  };

} // end namespace mu2e
#endif /* Mu2eG4_Mu2eWorld_hh */
