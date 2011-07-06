//
// Called at every G4 step.
//
// $Id: SteppingAction.cc,v 1.23 2011/07/06 22:45:34 logash Exp $
// $Author: logash $
// $Date: 2011/07/06 22:45:34 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <cstdio>
#include <cmath>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/SteppingAction.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
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

  SteppingAction::SteppingAction( const SimpleConfig& config ):

    // Default values for parameters that are optional in the run time configuration.
    _doKillLowEKine(false),
    _doKillInHallAir(false),
    _killerVerbose(false),
    _eKineMin(0.),
    _debugEventList(),
    _debugTrackList(),
    
    // Other parameters.
    _hallAirPhysVol(0),
    _lastPosition(),
    _lastMomentum(),
    _zref(0.),

    // Things related to time virtual detector
    _collection(0),
    _sizeLimit(config.getInt("g4.steppingActionStepsSizeLimit",0)),
    _currentSize(0),
    _simID(0),
    _productGetter(0) {

    // Look up parameter values in the run time configuration.
    _doKillLowEKine  = config.getBool("g4.killLowEKine",                _doKillLowEKine);
    _doKillInHallAir = config.getBool("g4SteppingAction.killInHallAir", _doKillInHallAir);
    _killerVerbose   = config.getBool("g4SteppingAction.killerVerbose", _killerVerbose);


    // If this cut is enabled, the cut value must be supplied in the run time config.
    // It is also used in StackingAction.
    if ( _doKillLowEKine ){
      _eKineMin = config.getDouble("g4.eKineMin");
    }

    vector<int> tmp1;
    config.getVectorInt( "g4.steppingActionEventList", tmp1, vector<int>() );
    _debugEventList.add(tmp1);

    vector<int> tmp2;
    config.getVectorInt( "g4.steppingActionTrackList", tmp2, vector<int>() );
    _debugTrackList.add(tmp2);

    //config.getVectorInt( "g4.steppingActionEventList", _debugEventList, vector<int>() );
    //config.getVectorInt( "g4.steppingActionTrackList", _debugTrackList, vector<int>() );
    // Get list of events for which to make debug printout.
    /*
    string key("g4.steppingActionEventList");
    if ( config.hasName(key) ){
      vector<int> list;
      config.getVectorInt(key,list);
      _debugEventList.add(list);
    }

    // Get list of tracks (within the above events) for which to make debug printout.
    // If empty list, then make printout for all tracks.
    string key2("g4.steppingActionTrackList");
    if ( config.hasName(key2) ){
      vector<int> list;
      config.getVectorInt(key2,list);
      _debugTrackList.add(list);
    }
    */

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
  void printit( G4String const& s,
                G4int id,
                G4ThreeVector const& pos,
                G4ThreeVector const& mom,
                double localTime,
                double globalTime ){

    // It is easier to line up printout in columns with printf than with cout.
    printf ( "%-8s %4d %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %13.4f %13.4f\n",
             s.data(), id,
             pos.x(), pos.y(), pos.z(),
             mom.x(), mom.y(), mom.z(),
             mom.mag(),
             localTime, globalTime);
  }

  void SteppingAction::beginRun(PhysicsProcessInfo& processInfo,
				CLHEP::Hep3Vector const& mu2eOrigin ){
    _hallAirPhysVol = getPhysicalVolumeOrThrow("HallAir");
    _currentSize    = 0;
    _processInfo    = &processInfo;
    _mu2eOrigin     =  mu2eOrigin;
  }

  void SteppingAction::UserSteppingAction(const G4Step* step){

    _nSteps++;

    G4Track* track = step->GetTrack();

    // Pre and post stepping points.
    G4StepPoint const* prept  = step->GetPreStepPoint();
    G4StepPoint const* postpt = step->GetPostStepPoint();

    // Save hits in time virtual detector
    for( unsigned int i=0; i<tvd_time.size(); ++i ) {
      if( prept->GetGlobalTime()<=tvd_time[i] && postpt->GetGlobalTime()>tvd_time[i] ) {
	addTimeVDHit(step,i+1);
      }
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

    if ( _doKillInHallAir &&  killInHallAir(track) ){
      killTrack( track, ProcessCode::mu2eHallAir, fStopAndKill);
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

    printf ( "\n");
    fflush(stdout);

  } // end UserSteppingAction

  // Kill tracks that drop below the kinetic energy cut.
  // It might be smarter to program G4 to do this itself?
  bool SteppingAction::killLowEKine ( const G4Track* trk ){
    if ( trk->GetKineticEnergy() < _eKineMin ){
      if ( _killerVerbose ){
        cout << "Killed track: low energy. " << trk->GetTrackID() << endl;
      }
      return true;
    }
    return false;
  }

  // Kill tracks that enter the hall air.
  bool SteppingAction::killInHallAir( const G4Track* trk ){

    // If we are not in the hall air, keep the track.
    if ( trk->GetVolume() != _hallAirPhysVol ) {
      return false;
    }

    if ( _killerVerbose ){
      cout << "Killed track: in Hall Air. " << trk->GetTrackID() << endl;
    }
    return true;
  }

  // Kill tracks that take too many steps.
  bool SteppingAction::killTooManySteps( const G4Track* track ){

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
  void SteppingAction::killTrack( G4Track* track, ProcessCode::enum_type code, G4TrackStatus status ){

    // Get user track informaton object from the track.
    G4VUserTrackInformation* info = track->GetUserInformation();
    UserTrackInformation* tinfo = (UserTrackInformation*)info;

    // Record why the track was killed.
    tinfo->setProcessCode(ProcessCode(code));

    // Kill the track
    track->SetTrackStatus(status);
  }

  void SteppingAction::BeginOfEvent(StepPointMCCollection& outputHits,
				    art::ProductID const& simID,
				    art::EDProductGetter const* productGetter ) {
    _nKilledStepLimit = 0;
    _collection  = &outputHits;
    _simID         = &simID;
    _productGetter = productGetter;
  }

  void SteppingAction::EndOfEvent() {
  }

  void SteppingAction::BeginOfTrack() {
    _nSteps = 0;
  }

  void SteppingAction::EndOfTrack() {
  }


  G4bool SteppingAction::addTimeVDHit(const G4Step* aStep, int id){

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in time virtual detector "
			     << _currentSize << endl;
      }
      return false;
    }

    // Which process caused this step to end?
    G4String const& pname  = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    ProcessCode endCode(_processInfo->findAndCount(pname));

    // The point's coordinates are saved in the mu2e coordinate system.
    _collection->
      push_back(StepPointMC(art::Ptr<SimParticle>( *_simID, aStep->GetTrack()->GetTrackID(), _productGetter ),
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

} // end namespace mu2e
