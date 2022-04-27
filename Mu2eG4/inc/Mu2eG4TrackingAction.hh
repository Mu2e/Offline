#ifndef Mu2eG4_Mu2eG4TrackingAction_hh
#define Mu2eG4_Mu2eG4TrackingAction_hh
//
// Steering routine for user tracking actions.
// If Mu2e needs many different user tracking actions, they
// should be called from this class.
//
// Original author Rob Kutschke
//

#include "CLHEP/Vector/ThreeVector.h"

#include "Geant4/G4UserTrackingAction.hh"

#include "Offline/MCDataProducts/inc/SimParticle.hh"

#include "Offline/Mu2eG4/inc/EventNumberList.hh"
#include "Offline/Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "cetlib/cpu_timer.h"

#include <map>
#include <string>

namespace mu2e {

  // Forward declarations in mu2e namespace
  class Mu2eG4SteppingAction;
  class Mu2eG4PerThreadStorage;
  class PhysicalVolumeHelper;
  class Mu2eG4ResourceLimits;
  namespace Mu2eG4Config { class Top; }


  class Mu2eG4TrackingAction: public G4UserTrackingAction{

  public:

    Mu2eG4TrackingAction(const Mu2eG4Config::Top& conf,
                   Mu2eG4SteppingAction *,
                   Mu2eG4PerThreadStorage *pts);

    // These methods are required by G4
    virtual void PreUserTrackingAction (const G4Track* trk);
    virtual void PostUserTrackingAction(const G4Track* trk);

    // All methods after here are Mu2e specific.

    // Do all things that need to be done at the beginning/end of an event.
    void beginEvent();

    void endEvent();

    // Record start and end points of each track created by G4.
    // Copy to the event data.
    void saveSimParticleStart(const G4Track* trk);
    void saveSimParticleEnd  (const G4Track* trk);

    // Receive persistent volume information and save it for the duration of the run.
    void beginRun( const PhysicalVolumeHelper* physVolHelper,
                   PhysicsProcessInfo* processInfo,
                   CLHEP::Hep3Vector const& mu2eOrigin );

    // Accessors for status information.
    unsigned        nG4Tracks() const { return _currentSize;}
    bool overflowSimParticles() const { return _overflowSimParticles; }
    unsigned nKilledByFieldPropagator() const { return _nKilledByFieldPropagator; }
    unsigned nKilledStepLimit() const { return _numKilledTracks; }


  private:

    typedef SimParticleCollection::key_type    key_type;
    typedef SimParticleCollection::mapped_type mapped_type;
    typedef std::map<key_type,mapped_type>     map_type;

    // Lists of events and tracks for which to enable debug printout.
    EventNumberList _debugList;

    // Non-owing pointer to the utility that translates between transient and persistent
    // representations of info about volumes.
    const PhysicalVolumeHelper* _physVolHelper;

    Mu2eG4PerThreadStorage *perThreadObjects_;

    // Origin of Mu2e Coordinate system in the G4 world system.
    CLHEP::Hep3Vector _mu2eOrigin;

    // Event timer.
    cet::cpu_timer _timer;

    // Information about SimParticles is collected in this map
    // during the operation of G4.  This is not persistent.
    map_type _transientMap;

    // Limit maximum size of the steps collection
    unsigned _sizeLimit;
    unsigned _currentSize;
    bool _overflowSimParticles;
    double _mcTrajectoryMomentumCut;
    double _saveTrajectoryMomentumCut;
    int    _mcTrajectoryMinSteps;
    unsigned _nKilledByFieldPropagator;
    unsigned _numKilledTracks;
    double _rangeToIgnore;

    const Mu2eG4ResourceLimits & _mu2elimits;

    // Non-owning pointer to stepping action; lifetime of pointee is one run.
    Mu2eG4SteppingAction * _steppingAction;

    // Non-owning pointer to the information about physical processes;
    // lifetime of pointee is one run.
    PhysicsProcessInfo *  _processInfo;
    bool _printTrackTiming;

    bool _stepLimitKillerVerbose;

    // Some helper functions.
    void insertOrThrow(std::pair<int,SimParticle> const& value);

    // If the track passes, the min hits cut and the momentum cut, add the
    // trajectory information to the output data product.
    void swapTrajectory( const G4Track* trk );

    // the muon specific decay proper time; it is ignored if set to a negative value
    double _muonPreAssignedDecayProperTime;
    // the maximum specific decay proper time; the specific time above excludes
    // the min max time use
    double _muonMinPreAssignedDecayProperTime;
    // the minimum specific decay proper time;
    double _muonMaxPreAssignedDecayProperTime;

  };

} // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4TrackingAction_hh */
