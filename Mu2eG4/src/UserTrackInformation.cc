//
// Mu2e specific information about one G4 track.
//
// $Id: UserTrackInformation.cc,v 1.3 2013/12/02 20:13:10 genser Exp $
// $Author: genser $
// $Date: 2013/12/02 20:13:10 $
//
// Original author Rob Kutschke
//

#include <iostream>

// Mu2e includes
#include "Mu2eG4/inc/UserTrackInformation.hh"

using namespace std;

namespace mu2e{

  UserTrackInformation::UserTrackInformation():
    G4VUserTrackInformation("Mu2eTrackInfo"),
    _forcedStop(false),
    _code(),
    _muCapCode(),
    _momDirection(),
    _kinEnergy(0)
  {}

  UserTrackInformation::~UserTrackInformation(){
  }

  void UserTrackInformation::Print() const{
  }

} // end namespace mu2e
