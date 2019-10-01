#ifndef Mu2eG4_Mu2eStudyWorld_hh
#define Mu2eG4_Mu2eStudyWorld_hh
//
// Construct the Mu2e G4 world and serve information about that
// world. Note that Mu2eStudyWorld inherits from Mu2eUniverse
//
// $Id: Mu2eStudyWorld.hh,v 1.1 2012/11/16 23:30:09 genser Exp $
// $Author: genser $
// $Date: 2012/11/16 23:30:09 $
//
// Original author KLG based on  Rob Kutschke Mu2eWorld
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
#include "Mu2eG4/inc/InitEnvToolBase.hh"
#include "Mu2eG4/inc/Mu2eUniverse.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "Mu2eG4/inc/FieldMgr.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "GeomPrimitives/inc/TubsParams.hh"

//G4 includes
#include "G4String.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"

namespace mu2e {

  // Forward references within mu2e namespace.
  class SensitiveDetectorHelper;

  class Mu2eStudyWorld : public Mu2eUniverse {
  public:

    Mu2eStudyWorld(const fhicl::ParameterSet& pset,
                   SensitiveDetectorHelper *sdHelper/*no ownership passing*/);

    Mu2eStudyWorld();
    ~Mu2eStudyWorld();

    // Construct everything.
    // The non-const return type is eventually required 
    // by G4VUserDetectorConstruction::Construct();
    //G4VPhysicalVolume * construct();
      
      virtual G4VPhysicalVolume * construct() override;
      
      virtual void constructSDandField() override;


  private:

    void constructStepLimiters();

    SensitiveDetectorHelper*          sdHelper_; // Non-owning
    std::unique_ptr<InitEnvToolBase>  constructEnv_;

    fhicl::ParameterSet               pset_;

    // _verbosityLevel in the base class
 
    bool        writeGDML_;
    std::string gdmlFileName_;
    std::string g4stepperName_;
    double      bfieldMaxStep_;
 

  };

} // end namespace mu2e
#endif /* Mu2eG4_Mu2eStudyWorld_hh */
