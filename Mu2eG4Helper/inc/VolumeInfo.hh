#ifndef Mu2eG4Helper_VolumeInfo_hh
#define Mu2eG4Helper_VolumeInfo_hh
//
// Information about a physical volume.  Used by Mu2eWorld and its utility routines.
// The center information is not fully general: it does not know about rotations
// and is useful only for the top few levels of the detector.
//
//
//
// Original author Rob Kutschke
//

#include <string>

#include "CLHEP/Vector/ThreeVector.h"

class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;

namespace mu2e {

  class Mu2eWorld;

  class VolumeInfo{

  public:

    // Depracated: this will go away.
    VolumeInfo():
      name(),
      solid(0),
      logical(0),
      physical(0),
      centerInParent(),
      centerInWorld(){}

    VolumeInfo( const std::string& pName);

    VolumeInfo( const std::string&       pName,
                const CLHEP::Hep3Vector& inParent,
                const CLHEP::Hep3Vector& parentInWorld);

    // Compiler written versions will be correct for:
    // destructor, copy constructor, assignment operator.

    // The name of this volume as known to G4.
    std::string name;

    // Non-owning pointers to volume information.
    G4VSolid*          solid;
    G4LogicalVolume*   logical;
    G4VPhysicalVolume* physical;

    // Location information in two coordinate systems.
    CLHEP::Hep3Vector      centerInParent;
    CLHEP::Hep3Vector      centerInWorld;

    CLHEP::Hep3Vector centerInMu2e() const { return centerInWorld - mu2eOriginInWorld(); }

  private:
    static const CLHEP::Hep3Vector& mu2eOriginInWorld();
  };

}

#endif /* Mu2eG4Helper_VolumeInfo_hh */
