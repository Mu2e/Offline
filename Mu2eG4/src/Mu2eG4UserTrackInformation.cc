//
// Mu2e specific information about one G4 track.
//
//
// Original author Rob Kutschke
//

#include <iostream>

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4UserTrackInformation.hh"

using namespace std;

namespace mu2e{

  Mu2eG4UserTrackInformation::Mu2eG4UserTrackInformation():
    G4VUserTrackInformation("Mu2eTrackInfo"),
    _forcedStop(false),
    _code(),
    _muCapCode(),
    _momDirection(),
    _kinEnergy(0)
  {}

  Mu2eG4UserTrackInformation::~Mu2eG4UserTrackInformation(){
  }

  void Mu2eG4UserTrackInformation::Print() const{
  }

} // end namespace mu2e
