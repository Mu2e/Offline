#ifndef MCDataProducts_StatusG4_hh
#define MCDataProducts_StatusG4_hh
//
// Status information about running G4 for one event.
//
//
// Original author Rob Kutschke
//
// Notes
// 1) If the _status data member is intended as a summary of the more detailed information available
//    in this class.  If it has a value of 0 then G4 completed with no signficant error or warning.
//    If it has a positive value then some issue arose during the run of G4 for this event;
//    the larger the number the worse the problem.  There are no negative values.
//    An enum has been used to define a large positive number that represents an undefined status;
//    only default constructed objects should have this status.
//
//    At this writing allowed non-zero values are:
//     1  - one or more tracks were stopped in SteppingAction because they took too many steps.
//          Usually these are particles trapped in the field. (or killed by the Field Propagator)
//     10 - The SimParticleCollection has overflowed.  So there can be hits that point back to
//          a Simparticle that is not in the collection.
//
// C++ includes
#include <iosfwd>
#include <climits>

namespace mu2e {

  class StatusG4 {

  public:

    // The undefined value is a large positive number.
    // Switch to numerical_limits<int>().max() in C++11
    enum enum_type { ok=0, undefined=INT_MAX };

    // This c'tor is required for ROOT.
    StatusG4():
      _status(undefined),
      _nG4Tracks(0),
      _overflowSimParticles(0),
      _nKilledStepLimit(0),
      _nKilledByFieldPropagator(0),
      _cpuTime(0),
      _realTime(0){}

    StatusG4( int status,
              int nG4Tracks,
              bool overflowSimParticles,
              int nKilledStepLimit,
              int nKilledByFieldPropagator,
              float cpuTime,
              float realTime
              ):
      _status(status),
      _nG4Tracks(nG4Tracks),
      _overflowSimParticles(overflowSimParticles),
      _nKilledStepLimit(nKilledStepLimit),
      _nKilledByFieldPropagator(nKilledByFieldPropagator),
      _cpuTime(cpuTime),
      _realTime(realTime)
    {}

    // Accept compiler written versions of d'tor, copy c'tor and assignment operator.

    void swap( StatusG4& );
    void print ( std::ostream&, bool newLine=true ) const;

    int                status() const { return _status; }
    int             nG4Tracks() const { return _nG4Tracks; }
    bool overflowSimParticles() const { return _overflowSimParticles; }
    int      nKilledStepLimit() const { return _nKilledStepLimit; }
    int nKilledByFieldPropagator() const { return _nKilledByFieldPropagator; }
    float             cpuTime() const { return _cpuTime; }
    float            realTime() const { return _realTime; }

    // Add the contents of another StatusG4 object to this one.
    // Used by event mixing to form a summary from multiple events.
    void add( StatusG4 const& );

  private:

    // Status=0 is all good.  Higher numbers indicate some issue has occured. See note 1.
    int _status;

    // Total number of G4 tracks that reach TrackingAction.
    // These are stored as SimParticles and that collection has a limited size; indicate an overflow.
    int  _nG4Tracks;
    bool _overflowSimParticles;

    // Number of G4Tracks killed in SteppingAction because there were too many steps.
    int _nKilledStepLimit;

    // Number of G4Tracks killed by the Field Propagator
    int _nKilledByFieldPropagator;

    // Execution time, in seconds.  CPU time is the sum of user time plus system time.
    // The times come from G4Timer - the least count is 10 ms and it appears that the resolution
    // is on the scale of the least count.
    float _cpuTime;
    float _realTime;

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   StatusG4 const& stat){
    stat.print(ost,false);
    return ost;
  }

} // end namespace mu2e

#endif /* MCDataProducts_StatusG4_hh */
