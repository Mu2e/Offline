#ifndef Mu2eUtilities_CosmicTrackUtils_hh
#define Mu2eUtilities_CosmicTrackUtils_hh
//
// Helper class to do work on CosmicTrack objects.
//

#include <tuple>

namespace mu2e {

  class CosmicTrack;

  // Function to make an std::tuple of kinkal track parameters from a CosmicTrack.
  std::tuple <double, double, double, double, double, double> KinKalTrackParams( CosmicTrack const& trk);

}

#endif
