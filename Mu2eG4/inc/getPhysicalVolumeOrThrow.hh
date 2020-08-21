#ifndef Mu2eG4_getPhysicalVolumeOrThrow_hh
#define Mu2eG4_getPhysicalVolumeOrThrow_hh
//
// Free function wrapper around:
//
//  G4PhysicalVolumeStore::GetInstance ()->GetVolume("HallAir");
//
// The wrapper does the job of throwing if the pointer comes back null.
//
//
// Original author Rob Kutschke
//
//

class G4VPhysicalVolume;
class G4String;

namespace mu2e {

  G4VPhysicalVolume* getPhysicalVolumeOrThrow( G4String const& name, bool mustHave=true );

} // end namespace mu2e
#endif /* Mu2eG4_getPhysicalVolumeOrThrow_hh */
