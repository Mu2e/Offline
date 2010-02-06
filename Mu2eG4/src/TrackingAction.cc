//
// Steering routine for user tracking actions. 
// If Mu2e needs many different user tracking actions, they
// should be called from this class.
//
// $Id: TrackingAction.cc,v 1.1 2010/02/06 19:39:09 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/02/06 19:39:09 $
//
// Original author Rob Kutschke
//

#include <iostream>

// Mu2e includes
#include "Mu2eG4/inc/TrackingAction.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

// G4 incldues
#include "globals.hh"
#include "G4RunManager.hh"

using namespace std;

namespace mu2e {

  TrackingAction::TrackingAction( const SimpleConfig& config):
    ncalls(0),
    debugList(){

    string name("g4.trackingActionEventList");
    if ( config.hasName(name) ){
      vector<int> list;
      config.getVectorInt(name,list);
      debugList.add(list);
    }

  }
  
  TrackingAction::~TrackingAction(){
  }

  void TrackingAction::PreUserTrackingAction(const G4Track* trk){
    if ( !debugList.inList() ) return;
    printInfo( trk, "Start new Track: ");

  }

  void TrackingAction::PostUserTrackingAction(const G4Track* trk){
    if ( !debugList.inList() ) return;
    printInfo( trk, "End Track:       ");

  }

  void TrackingAction::printInfo(const G4Track* trk, const string& text){

    const G4Event* event = G4RunManager::GetRunManager()->GetCurrentEvent();

    // Get some properties of the tracks.
    G4VPhysicalVolume* pvol = trk->GetVolume();
    G4String volName = (pvol !=0) ?
      pvol->GetName(): "Unknown Volume";

    G4ParticleDefinition* pdef = trk->GetDefinition();
    G4String partName = (pdef !=0) ?
      pdef->GetParticleName() : "Unknown Particle";

    cout << text
	 << setw(5) << event->GetEventID() << " "
	 << setw(4) << ncalls              << " "
	 << setw(4) << trk->GetTrackID()   << " "
	 << setw(4) << trk->GetParentID()  << " "
	 << setw(8) << partName            << " | "
	 << trk->GetPosition()             << " "
	 << trk->GetMomentum()             << " "
	 << volName
	 << endl;

  }



} // end namespace mu2e

