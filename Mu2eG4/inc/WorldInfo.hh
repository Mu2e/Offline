#ifndef WorldInfo_H
#define WorldInfo_H
//
//  A collection of information about the Mu2e G4 world that will
//  be interesting to others.
//
//
// $Id: WorldInfo.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include <string>

// Forward references.
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4Types.hh"
#include "G4ThreeVector.hh"


namespace mu2e {

  class WorldInfo {
  public:
    
    WorldInfo():
      worldSolid(0),
      worldLog(0),
      worldPhys(0),
      yEverest(0),
      Mu2eOrigin(),
      TrackerOrigin(){
    };
    ~WorldInfo(){}

    // Information about the world volume.
    G4Box*             worldSolid;
    G4LogicalVolume*   worldLog;
    G4VPhysicalVolume* worldPhys;

    // The top of the world.
    G4double           yEverest;

    // Location of the origin of the Mu2e coordinate system.
    G4ThreeVector      Mu2eOrigin;

    // Location of the center of the tracker.
    G4ThreeVector      TrackerOrigin;

  };

} // end namespace mu2e
#endif

