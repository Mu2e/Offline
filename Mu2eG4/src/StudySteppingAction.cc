//
// Called at every G4 step.
//
// $Id: StudySteppingAction.cc,v 1.8 2013/12/02 20:15:40 genser Exp $
// $Author: genser $
// $Date: 2013/12/02 20:15:40 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <cstdio>
#include <cmath>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Mu2eG4/inc/StudySteppingAction.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "Mu2eG4/inc/getPhysicalVolumeOrThrow.hh"

// G4 includes
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "G4String.hh"

using namespace std;

namespace mu2e {

  StudySteppingAction::StudySteppingAction( const SimpleConfig& config ):

    // Default values for parameters that are optional in the run time configuration.
    _doKillLowEKine(false),
    _killerVerbose(false),
    _eKineMin(0.),
    _killLowKineticEnergyPDG(),
    _eKineMinPDG(),
    _debugEventList(),
    _debugTrackList(),

    // Other parameters.
    _lastPosition(),
    _lastMomentum(),
    _zref(0.),
    _preStepEK(-1),

    // Things related to time virtual detector etc...

    _timeStepPointCollection(0),
    _timeStepPointCollectionSizeLimit(config.getInt("g4.steppingActionStepsSizeLimit",0)),
    _timeStepPointCollectionCurrentSize(0),

    _stepperStepPointCollection(0),
    _stepperStepPointCollectionSizeLimit(config.getInt("g4.steppingActionStepsSizeLimit",0)),
    _stepperStepPointCollectionCurrentSize(0),
    _saveSteps(config.getBool("g4.steppingActionSaveSteps", true)),
    // the above functionality is achieved in g4study2 by using
    // physics.producers.g4run.Mu2eG4SteppingOnlyCut : @local::stepperout

    _simID(0),
    _event(0) {

    // Look up parameter values in the run time configuration.
    _doKillLowEKine  = config.getBool("g4.killLowEKine",                _doKillLowEKine);
    config.getVectorString("g4SteppingAction.killInTheseVolumes", _killInTheseVolumes);
    _killerVerbose   = config.getBool("g4SteppingAction.killerVerbose", _killerVerbose);

    // this can be removed after a grace period
    // make sure the old job options are fixed:
    if(config.hasName("g4SteppingAction.killInHallAir")) {
      throw cet::exception("G4CONTROL")
        << "The parameter g4SteppingAction.killInHallAir "
        << "has been replaced with g4SteppingAction.killInTheseVolumes. "
        << "Please fix the geometry configuration file. "
        << "\n";
    }

    // If this cut is enabled, the cut value must be supplied in the run time config.
    // It is also used in StackingAction.
    if ( _doKillLowEKine ){
      _eKineMin = config.getDouble("g4.eKineMin");
      config.getVectorInt("g4.killLowEKinePDG", _killLowKineticEnergyPDG, vector<int>() );
      config.getVectorDouble("g4.eKineMinPDG", _eKineMinPDG, vector<double>() );
      if( _killLowKineticEnergyPDG.size() != _eKineMinPDG.size() ) {
        throw cet::exception("G4CONTROL")
          << "Sizes of g4.killLowEKinePDG and g4.eKineMinPDG do not match: "
          << _killLowKineticEnergyPDG.size() <<  " "
          << _eKineMinPDG.size() <<  " "
          << "\n";
      }
    }

    vector<int> tmp1;
    config.getVectorInt( "g4.steppingActionEventList", tmp1, vector<int>() );
    _debugEventList.add(tmp1);

    vector<int> tmp2;
    config.getVectorInt( "g4.steppingActionTrackList", tmp2, vector<int>() );
    _debugTrackList.add(tmp2);

    // Get maximum allowed number of steps per event
    _maxSteps = config.getInt("g4.steppingActionMaxSteps", 0);
    if( _maxSteps>0 ) {
      cout << "Limit maximum number of steps in SteppingAction to "
           << _maxSteps << endl;
    }

    // Get maximum global time
    _maxGlobalTime = config.getDouble("g4.steppingActionMaxGlobalTime", 0.0);
    if( _maxGlobalTime>0.1 ) {
      cout << "Limit maximum global time in SteppingAction to "
           << _maxGlobalTime << " ns" << endl;
    }

    // Read times for time virtual detector
    config.getVectorDouble( "g4.steppingActionTimeVD", tvd_time, vector<double>() );
    if( tvd_time.size()>0 ) {
      cout << "Time virtual detector is enabled. Particles are recorded at";
      for( unsigned int i=0; i<tvd_time.size(); ++i ) cout << " " << tvd_time[i];
      cout << " ns" << endl;
    }

  }

  // A helper function to manage the printout.
  void StudySteppingAction::printit( G4String const& s,
                                     G4int id,
                                     G4ThreeVector const& pos,
                                     G4ThreeVector const& mom,
                                     double localTime,
                                     double globalTime )
  {

    // It is easier to line up printout in columns with printf than with cout.
    printf ( "%-8s %4d %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %13.4f %13.4f\n",
             s.data(), id,
             pos.x(), pos.y(), pos.z(),
             mom.x(), mom.y(), mom.z(),
             mom.mag(),
             localTime, globalTime);
  }

  void StudySteppingAction::beginRun(PhysicsProcessInfo& processInfo,
                                CLHEP::Hep3Vector const& mu2eOrigin ){

    _timeStepPointCollectionCurrentSize    = 0;
    _stepperStepPointCollectionCurrentSize = 0;
    _processInfo    = &processInfo;
    _mu2eOrigin     =  mu2eOrigin;

  }

  void StudySteppingAction::finishConstruction(){

    for(unsigned i=0; i<_killInTheseVolumes.size(); ++i) {
      if ( true ){
        std::cout<<"Adding G4 killer volume = "<<_killInTheseVolumes[i]<<std::endl;
      }
      _killerVolumes.insert(getPhysicalVolumeOrThrow(_killInTheseVolumes[i]));
    }

  }

  void StudySteppingAction::UserSteppingAction(const G4Step* step){

    _nSteps++;

    G4Track* track = step->GetTrack();

    // Pre and post stepping points.
    G4StepPoint const* prept  = step->GetPreStepPoint();
    G4StepPoint const* postpt = step->GetPostStepPoint();

    // Get kinetic energy at the begin of the step
    _preStepEK = prept->GetKineticEnergy();

    G4VUserTrackInformation* info = track->GetUserInformation();
    UserTrackInformation* tinfo   = static_cast<UserTrackInformation*>(info);

    tinfo->setStepInfo(_preStepEK, _nSteps);

    // Save hits in time virtual detector
    for( unsigned int i=0; i<tvd_time.size(); ++i ) {
      if( prept->GetGlobalTime()<=tvd_time[i] && postpt->GetGlobalTime()>tvd_time[i] ) {
        addTimeVDHit(step,i+1);
      }
    }

    if(_saveSteps) {
      // Save the step point
      addStepPoint(step);
    }

    // Have we reached maximum allowed number of steps per track?
    if( _maxSteps>0 && _nSteps>_maxSteps ) {
      cout << "SteppingAction: kill particle pdg="
           << track->GetDefinition()->GetPDGEncoding()
           << " due to large number of steps." << endl;
      killTrack( track, ProcessCode::mu2eMaxSteps, fStopAndKill);
      ++_nKilledStepLimit;
    }

    // Have we reached maximum allowed global time?
    if( _maxGlobalTime>0.1 && track->GetGlobalTime()>_maxGlobalTime ) {
      cout << "SteppingAction: kill particle pdg="
           << track->GetDefinition()->GetPDGEncoding()
           << " because maximum global time is reached." << endl;
      killTrack( track, ProcessCode::mu2eMaxGlobalTime, fStopAndKill);
    }

    if ( killInTheseVolumes(track) ){
      killTrack( track, ProcessCode::mu2eKillerVolume, fStopAndKill);
    } else if ( _doKillLowEKine && killLowEKine(track) ){
      killTrack( track, ProcessCode::mu2eLowEKine, fStopAndKill);
    } else if ( _maxSteps>0 && killTooManySteps(track) ) {
      killTrack( track, ProcessCode::mu2eMaxSteps, fStopAndKill);
    }

    // Do we want to do make debug printout for this event?
    if ( !_debugEventList.inList() ) return;

    // Get information about this track.
    G4int id = track->GetTrackID();

    // If no tracks are listed, then printout for all tracks.
    // If some tracks are listed, then printout only for those tracks.
    if ( _debugTrackList.size() > 0 ){
      if ( !_debugTrackList.inList(id) ) return;
    }

    //G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();

    // Position and momentum at the the pre point.
    G4ThreeVector const& pos = prept->GetPosition();
    G4ThreeVector const& mom = prept->GetMomentum();

    // On the last step on a track the post step point does not have an
    // associated physical volume. So we need to protect against that.
    G4String preVolume, postVolume;
    G4int preCopy(-1), postCopy(-1);

    // Get the names if they are defined.
    if ( prept->GetPhysicalVolume() ){
      preVolume = prept->GetPhysicalVolume()->GetName();
      preCopy   = prept->GetPhysicalVolume()->GetCopyNo();
    }

    if ( postpt->GetPhysicalVolume() ){
      postVolume = postpt->GetPhysicalVolume()->GetName();
      postCopy   = postpt->GetPhysicalVolume()->GetCopyNo();
    }

    // On the forward trace, save the particle status at the start of the
    // last reporting volume.
    bool save = (std::abs(pos.z()+_zref)<0.0001) && (mom.z() > 0.);

    // On the backward trace, report the position when at the first
    // reporting volume.
    //bool report = (std::abs(pos.z()-_zref)<0.0001) && (mom.z() < 0.);

    // Save the status.
    if ( save ){
      _lastPosition = prept->GetPosition();
      _lastMomentum = prept->GetMomentum();
    }

    // Status report.
    printf ( "Step number: %d\n", _nSteps );

    printit ( "Pre: ", id,
              prept->GetPosition(),
              prept->GetMomentum(),
              prept->GetLocalTime(),
              prept->GetGlobalTime()
              );

    printit ( "Step:", id,
              track->GetPosition(),
              track->GetMomentum(),
              track->GetLocalTime(),
              track->GetGlobalTime()
              );

    printit ( "Post: ", id,
              postpt->GetPosition(),
              postpt->GetMomentum(),
              postpt->GetLocalTime(),
              postpt->GetGlobalTime()
              );
    fflush(stdout);

    cout << "Pre  Volume and copy: " << preVolume  << " " << preCopy  << endl;
    cout << "Post Volume and copy: " << postVolume << " " << postCopy << endl;

    G4String pname = Mu2eG4UserHelpers::findStepStoppingProcessName(step);
    cout << "Process which defined step: " << pname << " " << ProcessCode::findByName(pname) << endl; 

    printf ( "\n");
    fflush(stdout);

  } // end UserSteppingAction

  // Kill tracks that drop below the kinetic energy cut.
  // It might be smarter to program G4 to do this itself?
  bool StudySteppingAction::killLowEKine ( const G4Track* trk ){
    if ( trk->GetKineticEnergy() <= _eKineMin ){
      if ( _killerVerbose ){
        cout << "Killed track: low energy. " << trk->GetTrackID() << endl;
      }
      return true;
    }

    int pdg(trk->GetDefinition()->GetPDGEncoding());
    for( size_t i=0; i<_killLowKineticEnergyPDG.size(); ++i ) {
      if( _killLowKineticEnergyPDG[i] == pdg ) {
        if( trk->GetKineticEnergy() <= _eKineMinPDG[i] ){
          if ( _killerVerbose ){
            cout << "Killed track PDG " 
                 << pdg << " : low energy. " 
                 << trk->GetTrackID() << endl;
          }
          return true;
        }
        break;
      }
    }

    return false;
  }

  // Kill tracks that enter the selected volumes
  bool StudySteppingAction::killInTheseVolumes( const G4Track* trk ){

    KillerVolumesCache::const_iterator p = _killerVolumes.find(trk->GetVolume());
    if( p == _killerVolumes.end() ) {
      return false;
    }

    if ( _killerVerbose ){
      cout << "Killed track " << trk->GetTrackID()
           << " in volume " << (*p)->GetName()
           << endl;
    }
    return true;
  }

  // Kill tracks that take too many steps.
  bool StudySteppingAction::killTooManySteps( const G4Track* track ){

    if( _nSteps <= _maxSteps ) {
      return false;
    }

    if ( _killerVerbose ) {
      cout << "SteppingAction: kill particle pdg="
           << track->GetDefinition()->GetPDGEncoding()
           << " due to large number of steps." << endl;
    }
    ++_nKilledStepLimit;
    return true;
  }

  // Record why the track is to be killed, then kill it.
  void StudySteppingAction::killTrack( G4Track* track, ProcessCode::enum_type code, G4TrackStatus status ){

    // Get user track informaton object from the track.
    G4VUserTrackInformation* info = track->GetUserInformation();
    UserTrackInformation* tinfo   = static_cast<UserTrackInformation*>(info);

    // Record why the track was killed.
    tinfo->setProcessCode(ProcessCode(code));

    // Kill the track
    track->SetTrackStatus(status);
  }

  void StudySteppingAction::BeginOfEvent(StepPointMCCollection& timeOutputHits,
                                         StepPointMCCollection& stepperOutputStepPoints,
                                         art::ProductID const& simID,
                                         art::Event const& event ) {
    _nKilledStepLimit           = 0;
    _timeStepPointCollection    = &timeOutputHits;
    _stepperStepPointCollection = &stepperOutputStepPoints;
    _simID                      = &simID;
    _event                      = &event;
  }

  void StudySteppingAction::EndOfEvent() {
  }

  void StudySteppingAction::BeginOfTrack() {
    _nSteps = 0;
  }

  void StudySteppingAction::EndOfTrack() {
  }


  G4bool StudySteppingAction::addTimeVDHit(const G4Step* aStep, int id){

    ++_timeStepPointCollectionCurrentSize;

    if( _timeStepPointCollectionSizeLimit>0 && 
        _timeStepPointCollectionCurrentSize>_timeStepPointCollectionSizeLimit ) {
      if( (_timeStepPointCollectionCurrentSize - _timeStepPointCollectionSizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in time virtual detector "
                             << _timeStepPointCollectionCurrentSize << endl;
      }
      return false;
    }

    // Which process caused this step to end?
    G4String pname = Mu2eG4UserHelpers::findStepStoppingProcessName(aStep);
    ProcessCode endCode(_processInfo->findAndCount(pname));

    // The point's coordinates are saved in the mu2e coordinate system.
    _timeStepPointCollection->
      push_back(StepPointMC(art::Ptr<SimParticle>(*_simID, 
                                                  aStep->GetTrack()->GetTrackID(), 
                                                  _event->productGetter(*_simID) ),
                            id,
                            0,
                            0,
                            aStep->GetPostStepPoint()->GetGlobalTime(),
                            aStep->GetPostStepPoint()->GetProperTime(),
                            aStep->GetPostStepPoint()->GetPosition() - _mu2eOrigin,
                            aStep->GetPostStepPoint()->GetMomentum(),
                            aStep->GetStepLength(),
                            endCode
                            ));

    return true;

  }

  // we add step point using the pre step point info (like the sensitive detectors do)
  G4bool StudySteppingAction::addStepPoint(const G4Step* aStep){

    ++_stepperStepPointCollectionCurrentSize;

    if( _stepperStepPointCollectionSizeLimit>0 && 
        _stepperStepPointCollectionCurrentSize>_stepperStepPointCollectionSizeLimit ) {
      if( (_stepperStepPointCollectionCurrentSize - _stepperStepPointCollectionSizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of steps reached in  stepperStepPointCollection"
                             << _stepperStepPointCollectionCurrentSize << endl;
      }
      return false;
    }

    // Which process caused this step to end?
    G4String pname = Mu2eG4UserHelpers::findStepStoppingProcessName(aStep);
    ProcessCode endCode(_processInfo->findAndCount(pname));

    // The point's coordinates are saved in the mu2e coordinate system.
    _stepperStepPointCollection->
      push_back(StepPointMC(art::Ptr<SimParticle>(*_simID, 
                                                  aStep->GetTrack()->GetTrackID(), 
                                                  _event->productGetter(*_simID) ),
                            aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(),
                            aStep->GetTotalEnergyDeposit(),
                            aStep->GetNonIonizingEnergyDeposit(),
                            aStep->GetPreStepPoint()->GetGlobalTime(),
                            aStep->GetPreStepPoint()->GetProperTime(),
                            aStep->GetPreStepPoint()->GetPosition() - _mu2eOrigin,
                            aStep->GetPreStepPoint()->GetMomentum(),
                            aStep->GetStepLength(),
                            endCode
                            ));

    return true;

  }

} // end namespace mu2e
