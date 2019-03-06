#ifndef Mu2eG4_constructTracker_hh
#define Mu2eG4_constructTracker_hh
//
// Free functions to construct the TTracker.
//
// Original author Rob Kutschke
//


class G4LogicalVolume;

namespace mu2e{

  class VolumeInfo;
  class SimpleConfig;

  VolumeInfo constructTrackerv5( VolumeInfo const& parent,
                                 SimpleConfig const& config );

}  // end namespace mu2e

#endif /* Mu2eG4_constructTracker_hh */
