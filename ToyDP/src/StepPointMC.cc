//
// A persistable class representing a point that is on a track and
// is also inside, or on the boundary of, some G4 volume.  This can be
// used for saving points on the trajectory of the tracking and
// cosmic ray veto systems and for non-senstive material that we wish
// to record for purposes of debugging fitters.  We may need a different
// class to hold the corresponding information for calorimeters.
//
// $Id: StepPointMC.cc,v 1.7 2011/05/20 22:21:07 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/20 22:21:07 $
//
// Original author Rob Kutschke

// Mu2e includes
#include "ToyDP/inc/StepPointMC.hh"

using namespace std;

namespace mu2e {

  void StepPointMC::print( ostream& ost, bool doEndl ) const {

    ost << "  trackId: "                        << _trackId
        << "  volumeId: "                       << _volumeId
        << "  energy deposit: "                 << _totalEnergyDeposit
        << "  non ionizing energy deposit: "    << _nonIonizingEnergyDeposit
        << "  position: "                       << _position
        << "  momentum: "                       << _momentum
        << "  time: "                           << _time
        << "  proper time: "                    << _proper
        << "  step length: "                    << _stepLength
	<< "  end process: "                    << _endProcessCode;

    if ( doEndl ){
      ost << endl;
    }
  }

} // namespace mu2e
