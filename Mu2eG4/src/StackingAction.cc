#include "Mu2eG4/inc/StackingAction.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"

using namespace std;

namespace mu2e {

  StackingAction::StackingAction( const SimpleConfig& config ):
    ncalls(0),
    nevents(0),
    doCosmicKiller(false){

    // Get control info from run time configuration.
    doCosmicKiller = config.getBool("g4.doCosmicKiller",false);
  }

  StackingAction::~StackingAction(){
  }

  G4ClassificationOfNewTrack
  StackingAction::ClassifyNewTrack(const G4Track* trk){
    ++ncalls;

    // When we get a few different tests, we can decide how to
    // make a common interface.
    if ( doCosmicKiller ){
      if ( cosmicKiller(trk) ){
	return fKill;
      }
    }

    return fUrgent;
  }
  

  void StackingAction::NewStage(){
  }

  void StackingAction::PrepareNewEvent(){ 
    ncalls = 0;
    ++nevents;
  }

  bool StackingAction::cosmicKiller( const G4Track* trk){

    // Get some properties of the tracks.
    G4VPhysicalVolume* pvol = trk->GetVolume();
    G4String volName = (pvol !=0) ?
      pvol->GetName(): "Unknown Volume";
    
    G4ParticleDefinition* pdef = trk->GetDefinition();
    G4String partName = (pdef !=0) ?
      pdef->GetParticleName() : "Unknown Particle";

    const G4ThreeVector& ppos = trk->GetPosition();
    G4ThreeVector p3mom = trk->GetMomentum();

    // Magic numbers for illustrative purposes.
    // You can probably get much more aggressive than this.
    // Get ycut from geometry and pcut from run time config.
    static const double ycut = 800.;
    static const double pcut =  50.;

    // Decide if we want to kill the track.
    bool killit = ( ppos.y() > ycut && p3mom.mag() < pcut);

    // Printout about the decision.
    if ( nevents < 5 ) {
      G4String killString = killit ? "Kill it": "";
      cout << "Cosmic Killer: " 
	   << setw(4)  << nevents << " "
	   << setw(4)  << ncalls << " "
	   << setw(4)  << trk->GetTrackID() << " "
	   << setw(4)  << trk->GetParentID() <<  " "
	   << setw(8)  << partName << "   |   "
	   << ppos     << " "
	   << volName  << " "
	   << p3mom.mag() << "      "
	   << killString
	   << endl;
    }

    return killit;
  }

} // end namespace mu2e

