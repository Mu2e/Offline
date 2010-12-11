#ifndef ToyDP_StatusG4_hh
#define ToyDP_StatusG4_hh
// 
// Status information about running G4 for one event.
//
// $Id: StatusG4.hh,v 1.1 2010/12/11 00:31:03 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/12/11 00:31:03 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <iosfwd>

// Mu2e includes
#include "ToyDP/inc/GenId.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"

// Includes from external packages.
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"

namespace mu2e {

  class StatusG4 {

  public:

    // This c'tor is required for ROOT.
    StatusG4(){};

    StatusG4( int status,
              int nG4Tracks,
              bool overflowSimParticles,
              int nKilledStepLimit,
              float cpuTime,
              float realTime
              ):
      _status(status),
      _nG4Tracks(nG4Tracks),
      _overflowSimParticles(overflowSimParticles),
      _nKilledStepLimit(nKilledStepLimit),
      _cpuTime(cpuTime),
      _realTime(realTime)
    {}

    // Accept compiler written versions of d'tor, copy c'tor and assignment operator.

    void swap( StatusG4& );
    void print ( std::ostream& ) const;

    int                status() const { return _status; }
    int             nG4Tracks() const { return _nG4Tracks; }
    bool overflowSimParticles() const { return _overflowSimParticles; }
    int      nKilledStepLimit() const { return _nKilledStepLimit; }
    float             cpuTime() const { return _cpuTime; }
    float            realTime() const { return _realTime; }

    
  private:

    // Status=0 is all good.  Higher numbers indicate some issue has occured.
    // The user must check the other information to learn if the condition is a problem for his
    // work.
    int _status;

    // Total number of G4 tracks that reach TrackingAction.  
    // These are stored as SimParticles and that collection has a limited size; indicate an overflow.
    int  _nG4Tracks;
    bool _overflowSimParticles;

    // Number of G4Tracks killed in SteppingAction because there were too many steps.
    int _nKilledStepLimit;

    // Execution time, in seconds.  CPU time is the sum of user time plus system time.
    // The times come from G4Timer - not sure what is underneath that or the resolution.
    float _cpuTime;
    float _realTime;

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   StatusG4 const& stat){
    stat.print(ost);
    return ost;
  }


} // end namespace mu2e 

#endif
