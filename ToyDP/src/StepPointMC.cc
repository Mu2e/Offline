//
// A persistable class representing a point that is on a track and
// is also inside, or on the boundary of, some G4 volume.  This can be
// used for saving points on the trajectory of the tracking and 
// cosmic ray veto systems and for non-senstive material that we wish 
// to record for purposes of debugging fitters.  We may need a different 
// class to hold the corresponding information for calorimeters.
//
// $Id: StepPointMC.cc,v 1.2 2009/10/22 16:34:30 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/10/22 16:34:30 $
//
// Original author Rob Kutschke

// Mu2e incldues
#include "ToyDP/inc/StepPointMC.hh"

using namespace std;

namespace mu2e {

  void StepPointMC::print( ostream& ost, bool doEndl ) const {

    ost << "  trackId: "        << _trackId 
	<< "  volumeId: "       << _volumeId
	<< "  energy deposit: " << _edep
	<< "  position: "       << _position
	<< "  momentum: "       << _momentum
	<< "  time: "           << _time;

    if ( doEndl ){
      ost << endl;
    }
  }

} // namespace mu2e
