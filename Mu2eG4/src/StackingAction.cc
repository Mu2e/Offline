//
// Steering routine for user stacking actions. 
//
// $Id: StackingAction.cc,v 1.7 2011/02/13 22:33:10 logash Exp $
// $Author: logash $
// $Date: 2011/02/13 22:33:10 $
//
// Original author Rob Kutschke
//
// The actions steered by this class are:
//  1) cosmicKiller
//      - Kill low momentum secondaries of cosmic rays if those
//        secondaries are very unlikely to ever reach the detector. 
//      - The first implementation of this code is rather crude and
//        is presented to show techniques for accessing information.
//
//  Notes
//   1) Only one instance of this class should be instantiated in a job.
//   2) It is instantiated and registered with G4 in G4_plugin.cc
//   3) Once the G4 geometry has been created, the physical volumes
//      have static locations in memory.  So we can compare if two
//      physical volumes are the same object by comparing pointers,
//      rather than by comparing their names using string comparisons.
//   4) We get pointers to the physical volumes we care about in the
//      method PrepareNewEvent().  These pointer only change when the
//      G4 geometry changes but we have no way to detect that. So
//      do it every event.
//

// Mu2e includes
#include "Mu2eG4/inc/StackingAction.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "FWCore/Utilities/interface/Exception.h"
#include "Mu2eUtilities/inc/PDGCode.hh"

// G4 includes
#include "G4PhysicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"

using namespace std;

namespace mu2e {

  StackingAction::StackingAction( const SimpleConfig& config ):
    ncalls(0),
    nevents(0),
    doCosmicKiller(false),
    killLevel(0),
    _pdgToDrop(),
    dirtBodyPhysVol(0),
    dirtCapPhysVol(0){

    // Get control info from run time configuration.
    doCosmicKiller = config.getBool("g4.doCosmicKiller",false);
    killLevel = config.getInt("g4.cosmicKillLevel",0);

    // Get list of particles to keep or to drop in stepping action
    if ( config.hasName("g4.steppingActionDropPDG") ){
      config.getVectorInt("g4.steppingActionDropPDG",_pdgToDrop);
    }
    if( _pdgToDrop.size()>0 ) {
      cout << "Drop these particles in the SteppingAction: ";
      for( size_t i=0; i<_pdgToDrop.size(); ++i ) cout << _pdgToDrop[i] << ",";
      cout << endl;
    }

  }

  StackingAction::~StackingAction(){
  }

  G4ClassificationOfNewTrack
  StackingAction::ClassifyNewTrack(const G4Track* trk){
    ++ncalls;

    // When we get a few different tests, we can decide how to
    // make a more elegant way to call them all.
    if ( doCosmicKiller ){
      if ( cosmicKiller(trk) ){
        return fKill;
      }
    }

    if ( !_pdgToDrop.empty() ){
      if ( dropByPDGId(trk) ){
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

    // Find the addresses of some physical volumes of interest.
    // See notes 3 and 4.
    dirtBodyPhysVol = G4PhysicalVolumeStore::GetInstance ()->GetVolume("DirtBody");
    dirtCapPhysVol  = G4PhysicalVolumeStore::GetInstance ()->GetVolume("DirtCap");
  }

  bool StackingAction::cosmicKiller( const G4Track* trk){

    // Get some properties of the tracks.
    G4VPhysicalVolume* pvol = trk->GetVolume();
    G4ParticleDefinition* pdef = trk->GetDefinition();

    bool killit = false;

    if ( killLevel > 1 && pdef != 0 ) {
      int pType = pdef->GetPDGEncoding();

      // just kill anything electromagnetic in the dirt
      killit =  ( pvol == dirtBodyPhysVol &&
                  (pType == PDGCode::e_minus || 
                   pType ==  PDGCode::e_plus || 
                   PDGCode::gamma ) );
    } else {

      const G4ThreeVector& ppos = trk->GetPosition();
      G4ThreeVector p3mom = trk->GetMomentum();
      
      // Magic numbers for illustrative purposes.
      // You can probably get much more aggressive than this.
      // Get ycut from geometry and pcut from run time config.
      static const double ycut = 800.* CLHEP::mm;
      static const double pcut =  50.* CLHEP::MeV;
      
      // Decide if we want to kill the track.
      killit = ( ppos.y() > ycut && p3mom.mag() < pcut);
    }

    // Printout about the decision.
    if ( nevents < 20 ) {
      G4String volName = (pvol !=0) ? pvol->GetName(): "Unknown Volume";
      G4String partName = (pdef !=0) ?
        pdef->GetParticleName() : "Unknown Particle";
      G4String killString = killit ? "Kill it": "";
      const G4ThreeVector& ppos = trk->GetPosition();
      G4ThreeVector p3mom = trk->GetMomentum();

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

      // Check to see if we are in the dirt.
      if ( pvol == dirtBodyPhysVol ){
        cout << "Cosmic Killer: tag Dirt Body" << endl;
      } else if ( pvol == dirtCapPhysVol ) {
        cout << "Cosmic Killer: tag Dirt Cap" << endl;
      }

    }


    return killit;
  }

  // Return true if the particle Id of this track is in the list.
  bool StackingAction::dropByPDGId( G4Track const* track ){
    
    int id(track->GetDefinition()->GetPDGEncoding()); 
    for( size_t i=0; i<_pdgToDrop.size(); ++i ) {
      if( _pdgToDrop[i] == id ) {
        // cout << "Killing track from list: " << id << endl;
        return true;
      }
    }
    return false;
  }

} // end namespace mu2e

