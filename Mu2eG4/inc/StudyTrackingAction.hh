#ifndef Mu2eG4_StudyTrackingAction_hh
#define Mu2eG4_StudyTrackingAction_hh
//
// Steering routine for user tracking actions.
// If Mu2e needs many different user tracking actions, they
// should be called from this class.
//
// $Id: StudyTrackingAction.hh,v 1.5 2014/08/25 20:01:30 genser Exp $
// $Author: genser $
// $Date: 2014/08/25 20:01:30 $
//
// Original author Rob Kutschke
//

#include "CLHEP/Vector/ThreeVector.h"

#include "G4UserTrackingAction.hh"

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eG4/inc/EventNumberList.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "cetlib/cpu_timer.h"

#include <map>
#include <string>


namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;
  class StudySteppingAction;

  class StudyTrackingAction: public G4UserTrackingAction{

  public:

    StudyTrackingAction( const SimpleConfig& config, StudySteppingAction *);
    virtual ~StudyTrackingAction();

    // These methods are required by G4
    virtual void PreUserTrackingAction (const G4Track* trk);
    virtual void PostUserTrackingAction(const G4Track* trk);

    // All methods after here are Mu2e specific.

    // Do all things that need to be done at the beginning/end of an event.
    void beginEvent( art::Handle<GenParticleCollection> const& gensHandle, 
                     art::ProductID const& simID, 
                     art::Event const &event );
    void endEvent( SimParticleCollection& simParticles );

    // Record start and end points of each track created by G4.
    // Copy to the event data.
    void saveSimParticleStart(const G4Track* trk);
    void saveSimParticleEnd  (const G4Track* trk);

    // Receive persistent volume information and save it for the duration of the run.
    void beginRun( const PhysicalVolumeHelper& physVolHelper,
                   PhysicsProcessInfo& processInfo,
                   CLHEP::Hep3Vector const& mu2eOrigin );

    // Clean up at end of run.
    void endRun();

    // Accessors for status information.
    int             nG4Tracks() const { return _currentSize;}
    bool overflowSimParticles() const { return _overflowSimParticles; }

  private:

    typedef SimParticleCollection::key_type    key_type;
    typedef SimParticleCollection::mapped_type mapped_type;
    typedef std::map<key_type,mapped_type>     map_type;

    // Lists of events and tracks for which to enable debug printout.
    EventNumberList _debugList;

    // Non-owing pointer to the utility that translates between transient and persistent
    // representations of info about volumes.
    const PhysicalVolumeHelper* _physVolHelper;

    // Origin of Mu2e Coordinate system in the G4 world system.
    CLHEP::Hep3Vector _mu2eOrigin;

    // Event timer.
    cet::cpu_timer _timer;

    // Information about SimParticles is collected in this map
    // during the operation of G4.  This is not persistent.
    map_type _transientMap;

    // Limit maximum size of the steps collection
    int _sizeLimit;
    int _currentSize;
    bool _overflowSimParticles;
    double _saveTrajectoryMomentumCut;

    // Non-owning pointer to stepping action; lifetime of pointee is one run.
    StudySteppingAction * _steppingAction;

    // Non-owning pointer to the information about physical processes;
    // lifetime of pointee is one run.
    PhysicsProcessInfo *  _processInfo;
    bool _printTrackTiming;

    // Handle to the GenParticle collection; needed to make Ptrs into that collection.
    art::Handle<GenParticleCollection> const * _gensHandle;

    // Product ID of the SimParticleCollection; needed to make PTrs into that collection.
    art::ProductID _simID;

    // Non-owning pointer to the event; needed for some Ptr stuff.
    art::Event const * _event;

    // Some helper functions.
    void insertOrThrow(std::pair<int,SimParticle> const& value);

  };

} // end namespace mu2e

#endif /* Mu2eG4_StudyTrackingAction_hh */

