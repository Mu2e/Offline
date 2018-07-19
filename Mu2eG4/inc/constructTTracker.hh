#ifndef Mu2eG4_constructTTracker_hh
#define Mu2eG4_constructTTracker_hh
//
// Free functions to construct the TTracker.
//
// Original author Rob Kutschke
//


class G4LogicalVolume;

namespace mu2e{

  class VolumeInfo;
  class SimpleConfig;

  VolumeInfo constructTTrackerv5( VolumeInfo const& parent,
                                 SimpleConfig const& config );

}  // end namespace mu2e

#endif /* Mu2eG4_constructTTracker_hh */
