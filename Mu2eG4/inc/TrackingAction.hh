#ifndef Mu2eG4_TrackingAction_hh
#define Mu2eG4_TrackingAction_hh
//
// Steering routine for user tracking actions.
// If Mu2e needs many different user tracking actions, they
// should be called from this class.
//
// $Id: TrackingAction.hh,v 1.14 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <string>
#include <map>

// Framework includes
#include "art/Utilities/CPUTimer.h"

// Mu2e includes
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/EventNumberList.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"

// G4 includes.
#include "G4UserTrackingAction.hh"

// Other includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;
  class SteppingAction;

  class TrackingAction: public G4UserTrackingAction{

  public:

    TrackingAction( const SimpleConfig& config, SteppingAction *);
    virtual ~TrackingAction();

    // These methods are required by G4
    virtual void PreUserTrackingAction (const G4Track* trk);
    virtual void PostUserTrackingAction(const G4Track* trk);

    // All methods after here are Mu2e specific.

    // Do all things that need to be done at the beginning/end of an event.
    void beginEvent();
    void endEvent( SimParticleCollection& simParticles );

    // Record start and end points of each track created by G4.
    // Copy to the event data.
    void saveSimParticleStart(const G4Track* trk);
    void saveSimParticleEnd  (const G4Track* trk);

    // Receive persistent volume information and save it for the duration of the run.
    void beginRun( const PhysicalVolumeHelper& physVolHelper,
                   CLHEP::Hep3Vector const& mu2eOrigin );

    // Clean up at end of run.
    void endRun();

    // Check consistency of mother-daughter pointers.
    bool checkCrossReferences( bool doPrint, bool doThrow);

    // Accessors for status information.
    int             nG4Tracks() const { return _currentSize;}
    bool overflowSimParticles() const { return _overflowSimParticles; }


  private:

    typedef SimParticleCollection::key_type key_type;
    typedef SimParticleCollection::map_type map_type;

    // Lists of events and tracks for which to enable debug printout.
    EventNumberList _debugList;

    // Utility to translate between transient and persistent representations.
    const PhysicalVolumeHelper* _physVolHelper;

    CLHEP::Hep3Vector _mu2eOrigin;

    art::CPUTimer _timer;

    // Information about SimParticles is collected in this map
    // during the operation of G4.  This is not persistent.
    map_type _transientMap;

    // Debug printout.
    void printInfo(const G4Track* trk, const std::string& text, bool isEnd=false);

    // Limit maximum size of the steps collection
    int _sizeLimit;
    int _currentSize;
    bool _overflowSimParticles;

    // Non-owning pointer to stepping action.
    SteppingAction * _steppingAction;

    // Non-owning pointer to the information about physical processes.
    PhysicsProcessInfo  _processInfo;

    // Do we print the summary of the process information at the end of the job.
    bool _printPhysicsProcessSummary;

    // Control the saving of trajectories.
    // The first method does the big picture bookkeeping.
    // The second method decides yes/no for storing the trajectory of one track.
    void controlTrajectorySaving( const G4Track* trk);
    bool saveThisTrajectory( const G4Track* trk );

    // Some helper functions.
    void insertOrThrow(std::pair<int,SimParticle> const& value);
    G4String findStoppingProcess(G4Track const* track);
    ProcessCode findCreationCode(G4Track const* track);

  };

} // end namespace mu2e

#endif /* Mu2eG4_TrackingAction_hh */

