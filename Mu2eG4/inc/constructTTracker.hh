#ifndef constructTTracker_HH
#define constructTTracker_HH
//
// Free functions to construct various versions of the TTracker.
//
// $Id: constructTTracker.hh,v 1.1 2010/04/18 00:33:02 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/18 00:33:02 $
//
// Original author Rob Kutschke
//

#include "Mu2eG4/inc/VolumeInfo.hh"

class G4LogicalVolume;

namespace mu2e{

  class SimpleConfig;

  VolumeInfo constructTTrackerv1( G4LogicalVolume* mother,
                                  double zOff,
                                  SimpleConfig const& config );


}  // end namespace mu2e

#endif
