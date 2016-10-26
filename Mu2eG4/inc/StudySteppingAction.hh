#ifndef Mu2eG4_StudySteppingAction_hh
#define Mu2eG4_StudySteppingAction_hh
//
// Called at every G4 step.
//
// $Id: StudySteppingAction.hh,v 1.2 2012/12/17 16:59:30 gandr Exp $
// $Author: gandr $
// $Date: 2012/12/17 16:59:30 $
//
// Original author Rob Kutschke
//
#include <vector>
#include <string>

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
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

// Forward declarations outside of mu2e namespace.
class G4VPhysicalVolume;
class G4Track;

namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;

  class StudySteppingAction : public G4UserSteppingAction{

  public:
    StudySteppingAction( const SimpleConfig& config );
    ~StudySteppingAction(){}

    void UserSteppingAction(const G4Step*);

    void BeginOfEvent(StepPointMCCollection& timeOutputHits,
                      StepPointMCCollection& stepperOutputStepPoints,
                      art::ProductID const& simID,
                      art::Event const &event );
    void EndOfEvent();

    void BeginOfTrack();
    void EndOfTrack();

    int nKilledStepLimit() const { return _nKilledStepLimit; }

    // Called by G4_plugin.
    void beginRun(PhysicsProcessInfo&, CLHEP::Hep3Vector const& mu2eOrigin );

    // Called by G4_plugin: the final phase of the c'tor cannot be completed until after
    // G4 has initialized itself.
    void finishConstruction();

    G4ThreeVector const& lastPosition() const { return _lastPosition; }
    G4ThreeVector const& lastMomentum() const { return _lastMomentum; }

    void setZRef( G4double zref){
      _zref=zref;
    }
    // A helper function to manage the printout.
    static void printit( G4String const& s,
                         G4int id,
                         G4ThreeVector const& pos,
                         G4ThreeVector const& mom,
                         double localTime,
                         double globalTime );

  private:

    // Start: information from the run time configuration.

    // Which killers will be enabled?
    bool _doKillLowEKine;
    // A potentially emtpy list of volumes names where tracks are killed
    std::vector<std::string> _killInTheseVolumes;
    bool _killerVerbose;

    // Minimum energy cut.
    double _eKineMin;
    std::vector<int> _killLowKineticEnergyPDG;
    std::vector<double> _eKineMinPDG;

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

    // cached pointers to physical volumes on the _killInTheseVolumes list
    typedef std::set<const G4VPhysicalVolume*> KillerVolumesCache;
    KillerVolumesCache _killerVolumes;

    // Used in the turn around code.
    G4ThreeVector _lastPosition;
    G4ThreeVector _lastMomentum;
    G4double      _zref;

    // Kinetic energy at the beginning of the step point
    double _preStepEK;

    // Functions to decide whether or not to kill tracks.
    bool killTooManySteps ( const G4Track* );
    bool killLowEKine     ( const G4Track* );
    bool killInTheseVolumes( const G4Track* );

    // A helper function to kill the track and record the reason for killing it.
    void killTrack( G4Track* track, ProcessCode::enum_type code, G4TrackStatus status );

    // Collection to hold time virtual detector hits
    StepPointMCCollection* _timeStepPointCollection;

    // Limit maximum size of the steps collections
    int _timeStepPointCollectionSizeLimit;
    int _timeStepPointCollectionCurrentSize;

    // Collection to hold stepper points
    StepPointMCCollection* _stepperStepPointCollection;

    // Limit maximum size of the steps collections
    int _stepperStepPointCollectionSizeLimit;
    int _stepperStepPointCollectionCurrentSize;
    bool _saveSteps;

    // Information about the SimParticleCollection, needed to instantiate art::Ptr.
    art::ProductID const*  _simID;
    art::Event     const*  _event;

    // Non-owning pointer to the information about physical processes;
    // lifetime of pointee is one run.
    PhysicsProcessInfo *  _processInfo;

    // Origin of Mu2e Coordinate system in the G4 world system.
    CLHEP::Hep3Vector _mu2eOrigin;

    // List of times for time virtual detector
    std::vector<double> tvd_time;

    // Add time virtual detector hit to the collection
    G4bool addTimeVDHit(const G4Step*, int);

    // Add step point to the collection
    G4bool addStepPoint(const G4Step*);

  };

} // end namespace mu2e
#endif /* Mu2eG4_StudySteppingAction_hh */
