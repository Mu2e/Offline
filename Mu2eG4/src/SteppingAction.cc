//
// Called at every G4 step.
//
// $Id: SteppingAction.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include <cstdio>
#include <cmath>

#include "Mu2eG4/inc/SteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4String.hh"

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


void SteppingAction::UserSteppingAction(const G4Step* step)
{  
  /*
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

  if ( save   ) G4cout << "Save point: " << G4endl;
  if ( report ) G4cout << "Report point: " << G4endl;

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

  printf ( "\n");
  */

}

}

