//
// Called at every G4 step.
//
// $Id: SteppingAction.cc,v 1.7 2010/08/30 22:21:13 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/08/30 22:21:13 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <cstdio>
#include <cmath>

// Mu2e includes
#include "Mu2eG4/inc/SteppingAction.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4String.hh"


using namespace std;

namespace mu2e {

  SteppingAction::SteppingAction( const SimpleConfig& config ):
    _lastPosition(),
    _lastMomentum(),
    _zref(),
    _debugEventList(),
    _debugTrackList(){ 

    // Get list of events for which to make debug printout.
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

  }

  // A helper function to manage the printout.
  void printit( G4String const& s, 
                G4int id,
                G4ThreeVector const& pos,
                G4ThreeVector const& mom){

    // It is easier to line up printout in columns with printf than with cout.
    printf ( "%-8s %4d %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f\n",
             s.data(), id,
             pos.x(), pos.y(), pos.z(), 
             mom.x(), mom.y(), mom.z(),
             mom.mag());
  }


  void SteppingAction::UserSteppingAction(const G4Step* step){  
    
    // Do we want to do make debug printout for this event?
    if ( !_debugEventList.inList() ) return;

    // Get information about this track.
    G4Track* track = step->GetTrack();
    G4int id       = track->GetTrackID();

    // If no tracks are listed, then printout for all tracks.
    // If some tracks are listed, then printout only for those tracks.
    if ( _debugTrackList.size() > 0 ){
      if ( !_debugTrackList.inList(id) ) return;
    }

    G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();
    int eventNo = event->GetEventID();

    // Pre and post stepping points.
    G4StepPoint const* prept  = step->GetPreStepPoint();
    G4StepPoint const* postpt = step->GetPostStepPoint();

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
    bool report = (std::abs(pos.z()-_zref)<0.0001) && (mom.z() < 0.);

    // Save the status.
    if ( save ){
      _lastPosition = prept->GetPosition();
      _lastMomentum = prept->GetMomentum();
    }

    //  if ( save   ) G4cout << "Save point: " << G4endl;
    //if ( report ) G4cout << "Report point: " << G4endl;

    // Status report.
    printit ( "Pre: ", id, 
              prept->GetPosition(),
              prept->GetMomentum()
              );


    printit ( "Step:", id, 
              track->GetPosition(),
              track->GetMomentum()
              );


    printit ( "Post: ", id, 
              postpt->GetPosition(),
              postpt->GetMomentum()
              );
    fflush(stdout);

    cout << "Pre  Volume and copy: " << preVolume  << " " << preCopy  << endl;
    cout << "Post Volume and copy: " << postVolume << " " << postCopy << endl;

    printf ( "\n");
    fflush(stdout);

  }

} // end namespace mu2e

