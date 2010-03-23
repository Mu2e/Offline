//
// Called at every G4 step.
//
// $Id: SteppingAction.cc,v 1.3 2010/03/23 20:43:14 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/03/23 20:43:14 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <cstdio>
#include <cmath>

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"
#include "Mu2eG4/inc/SteppingAction.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4String.hh"


using namespace std;

namespace mu2e {

  SteppingAction::SteppingAction():
    _zref(){ 
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

    G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();
    int eventNo = event->GetEventID();


    // Build list of interesting events.
    /*
      static int const nadded(12);
      static int const nmissing(4);
      static int const added[nadded]     = {  0, 15, 34,  49, 61, 66, 74, 99, 128, 142, 164, 172};
      static int const missing[nmissing] = { 25, 41, 63, 144 };
      static EventNumberList  add( nadded,   added  );
      static EventNumberList miss( nmissing, missing);
    */

    static EventNumberList  add;
    static EventNumberList miss;

    // Skip uninteresting events.
    bool inList = ( add.inList() || miss.inList() );
    if ( !inList ) return;

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
    G4Track* track = step->GetTrack();

    // Status report.
    G4int id = step->GetTrack()->GetTrackID();
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

    cout << "Pre  Volume and copy: " << preVolume  << " " << preCopy  << endl;
    cout << "Post Volume and copy: " << postVolume << " " << postCopy << endl;

    printf ( "\n");

  }

} // end namespace mu2e

