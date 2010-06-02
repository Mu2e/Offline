//
// Called at every G4 step.
//
// $Id: SteppingAction.cc,v 1.5 2010/06/02 04:01:53 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/06/02 04:01:53 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <cstdio>
#include <cmath>
#include <set>
#include <string>

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
    _zref(){ 

    // Get list of events for which to make debug printout.
    string key("g4.steppingActionEventList");
    if ( config.hasName(key) ){
      vector<int> list;
      config.getVectorInt(key,list);
      _debugList.add(list);
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

    static set<G4VPhysicalVolume*> vnames;
    G4StepPoint const* prept2 = step->GetPreStepPoint();
    G4VPhysicalVolume* prevol2 = prept2->GetPhysicalVolume();
    if ( prevol2 ){
      G4String preVol = prept2->GetPhysicalVolume()->GetName();
      int preCopy     = prept2->GetPhysicalVolume()->GetCopyNo();
      string spreVol(preVol);
        if ( spreVol.find("Straw") == string::npos ) {
          if ( vnames.find(prevol2) == vnames.end() ){
            vnames.insert(prevol2);
            G4AffineTransform const& toLocal = step->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform();
            G4AffineTransform        toWorld = toLocal.Inverse();
            G4ThreeVector localOrigin(0.,0.,0.);
            G4ThreeVector worldOrigin = toWorld.TransformPoint(localOrigin);
            cout << "Local Origin for "
                 << setw(40) << preVol << ", "
                 << setw(4) << preCopy << " : "
                 << worldOrigin
                 << endl;
            /*
            if ( spreVol == "ProductionTarget" ){
              G4ThreeVector localxhat(1.,0.,0.);
              G4ThreeVector localyhat(0.,1.,0.);
              G4ThreeVector localzhat(0.,0.,1.);
              G4ThreeVector worldxhat = toWorld.TransformPoint(localxhat) - worldOrigin;
              G4ThreeVector worldyhat = toWorld.TransformPoint(localyhat) - worldOrigin;
              G4ThreeVector worldzhat = toWorld.TransformPoint(localzhat) - worldOrigin;
              cout << "ProtonTarget xhat: " << worldxhat << " " << worldxhat.mag() << endl;
              cout << "ProtonTarget yhat: " << worldyhat << " " << worldyhat.mag() << endl;
              cout << "ProtonTarget zhat: " << worldzhat << " " << worldzhat.mag() << endl;
            }
            */
          }
        }
    }

    // Do we want to do make debug printout for this event?
    if ( !_debugList.inList() ) return;

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
    fflush(stdout);

    cout << "Pre  Volume and copy: " << preVolume  << " " << preCopy  << endl;
    cout << "Post Volume and copy: " << postVolume << " " << postCopy << endl;

    printf ( "\n");
    fflush(stdout);

  }

} // end namespace mu2e

