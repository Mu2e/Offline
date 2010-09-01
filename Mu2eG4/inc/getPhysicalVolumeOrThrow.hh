#ifndef getPhysicalVolumeOrThrow_H
#define getPhysicalVolumeOrThrow_H 1
//
// Free function wrapper around:
//
//  G4PhysicalVolumeStore::GetInstance ()->GetVolume("HallAir");
//
// The wrapper does the job of throwing if the pointer comes back null.
//
// $Id: getPhysicalVolumeOrThrow.hh,v 1.1 2010/09/01 18:57:19 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/09/01 18:57:19 $
//
// Original author Rob Kutschke
//
//

class G4VPhysicalVolume;
class G4String;

namespace mu2e {
  
  G4VPhysicalVolume* getPhysicalVolumeOrThrow( G4String const& name, bool mustHave=true );

} // end namespace mu2e
#endif
