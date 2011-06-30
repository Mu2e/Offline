#ifndef Mu2eG4_SteppingAction_hh
#define Mu2eG4_SteppingAction_hh
//
// Called at every G4 step.
//
// $Id: SteppingAction.hh,v 1.14 2011/06/30 20:27:53 logash Exp $
// $Author: logash $
// $Date: 2011/06/30 20:27:53 $
//
// Original author Rob Kutschke
//
#include <vector>

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"
#include "Mu2eG4/inc/UserTrackInformation.hh"
#include "MCDataProducts/inc/ProcessCode.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

// G4 includes
#include "CLHEP/Vector/ThreeVector.h"
#include "G4UserSteppingAction.hh"
#include "G4TrackStatus.hh"
#include "G4ThreeVector.hh"

// Art includes
#include "art/Persistency/Provenance/ProductID.h"
#include "art/Persistency/Common/EDProductGetter.h"

// Forward declarations outside of mu2e namespace.
class G4VPhysicalVolume;
class G4Track;

namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;

  class SteppingAction : public G4UserSteppingAction{

  public:
    SteppingAction( const SimpleConfig& config );
    ~SteppingAction(){}

    void UserSteppingAction(const G4Step*);

    void BeginOfEvent(StepPointMCCollection& outputHits,
		      art::ProductID const& simID,
		      art::EDProductGetter const* productGetter );
    void EndOfEvent();

    void BeginOfTrack();
    void EndOfTrack();

    int nKilledStepLimit() const { return _nKilledStepLimit; }

    // Called by G4_plugin.
    void beginRun(PhysicsProcessInfo&, CLHEP::Hep3Vector const& mu2eOrigin );

    G4ThreeVector const& lastPosition() const { return _lastPosition; }
    G4ThreeVector const& lastMomentum() const { return _lastMomentum; }

    void setZRef( G4double zref){
      _zref=zref;
    }
    
  private:

    // Start: information from the run time configuration.

    // Which killers will be enabled?
    bool _doKillLowEKine;
    bool _doKillInHallAir;
    bool _killerVerbose;

    // Minimum energy cut.
    double _eKineMin;

    // Maximum allowed number of steps per event
    int _maxSteps;
    int _nSteps;
    int _nKilledStepLimit;

    // Maximum global time of particle
    double _maxGlobalTime;

    // Lists of events and tracks for which to enable debug printout.
    EventNumberList _debugEventList;
    EventNumberList _debugTrackList;

    // End: information from the run time configuration.

    // Address of the physical volume object for the hall air.
    G4VPhysicalVolume* _hallAirPhysVol;

    // Used in the turn around code.
    G4ThreeVector _lastPosition;
    G4ThreeVector _lastMomentum;
    G4double      _zref;

    // Functions to decide whether or not to kill tracks.
    bool killTooManySteps ( const G4Track* );
    bool killLowEKine     ( const G4Track* );
    bool killInHallAir    ( const G4Track* );

    // A helper function to kill the track and record the reason for killing it.
    void killTrack( G4Track* track, ProcessCode::enum_type code, G4TrackStatus status );

    // Collection to hold time virtual detector hits
    StepPointMCCollection* _collection;

    // Limit maximum size of the steps collection
    int _sizeLimit;
    int _currentSize;

    // Information about the SimParticleCollection, needed to instantiate art::Ptr.
    art::ProductID const *      _simID;
    art::EDProductGetter const* _productGetter;

    // Non-owning pointer to the information about physical processes;
    // lifetime of pointee is one run.
    PhysicsProcessInfo *  _processInfo;

    // Origin of Mu2e Coordinate system in the G4 world system.
    CLHEP::Hep3Vector _mu2eOrigin;

    // List of times for time virtual detector
    std::vector<double> tvd_time;

    // Add time virtual detector hit to the collection
    G4bool addTimeVDHit(const G4Step*, int);

  };

} // end namespace mu2e
#endif /* Mu2eG4_SteppingAction_hh */
