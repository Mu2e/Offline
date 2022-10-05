#ifndef Mu2eUtilities_TrackerBFieldInfo_hh
#define Mu2eUtilities_TrackerBFieldInfo_hh
//
// Find the magnetic field map that contains the tracker origin.
// Return the map on request.
//
// I made this a class so that we can add other accessors that
// transform the information received from the map.
//
// If a map cannot be located it throws an exception.
// If there is more than one map among the BField innerMaps, it
// also throws and exception.
//
// To use this to look at limits of a map:
// TrackerBFieldInfo info
//  cout << info.map().xmin() << " " << info.map().xmax() << endl;  // and so on
//  cout << info.map().getKey() << endl;                            // name of the map that was selected
//
// Original author Rob Kutschke
//

#include "Offline/BFieldGeom/inc/BFMap.hh"

#include <iosfwd>

namespace mu2e {


  class TrackerBFieldInfo {
  public:

    TrackerBFieldInfo();

    BFMap const& map() const { return *_map; }

  private:

    BFMap const* _map = nullptr;

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_TrackerBFieldInfo_hh */
