//
// Steering routine for user stacking actions.
//
// $Id: StackingAction.cc,v 1.17 2011/05/18 14:21:44 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/18 14:21:44 $
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
//   4) Tracks created in PrimaryGeneratorAction do not have a defined
//      volume pointer at stacking time.  Tracks created by G4 processes
//      do a have defined volume pointer at stacking time.

// Mu2e includes
#include "Mu2eG4/inc/StackingAction.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "cetlib/exception.h"
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
    _ncalls(0),
    _nevents(0),
    _doCosmicKiller(false),
    _killLevel(0),
    _cosmicpcut(0),
    _yaboveDirtYmin(0),
    _primaryOnly(false),
    _killLowKineticEnergy(false),
    _eKineMin(0.),
    _pdgToDrop(),
    _pdgToKeep(),
    _dirtBodyPhysVol(0),
    _dirtCapPhysVol(0),
    _dirtG4Ymin(0),
    _dirtG4Ymax(0){

    // Get control info from run time configuration.
    _doCosmicKiller       = config.getBool  ("g4.doCosmicKiller",   _doCosmicKiller );
    _killLevel            = config.getInt   ("g4.cosmicKillLevel",  _killLevel      );
    _cosmicpcut           = config.getDouble("g4.cosmicPcut",       _cosmicpcut     );
    _yaboveDirtYmin       = config.getDouble("g4.yaboveDirtYmin",   _yaboveDirtYmin );
    _primaryOnly          = config.getBool  ("g4.stackPrimaryOnly", _primaryOnly    );
    _killLowKineticEnergy = config.getBool  ("g4.killLowEKine",     _killLowKineticEnergy );

    config.getVectorInt("g4.stackingActionDropPDG", _pdgToDrop, vector<int>() );
    config.getVectorInt("g4.stackingActionKeepPDG", _pdgToKeep, vector<int>() );

    // If this cut is enabled, the cut value must be supplied in the run time config.
    // It is also used in SteppingAction.
    if ( _killLowKineticEnergy ){
      _eKineMin = config.getDouble("g4.eKineMin");
    }

    if ( !_pdgToDrop.empty() && !_pdgToKeep.empty() ){
      throw cet::exception("G4CONTROL")
        << "Both g4.stackingActionKeepPDG and g4.stackingActionDropPDG have entries: "
        << _pdgToDrop.size() <<  " "
        << _pdgToKeep.size() <<  " "
        << "\n";
    }

    if( _pdgToDrop.size()>0 ) {
      cout << "Drop these particles in the StackingAction: ";
      for( size_t i=0; i<_pdgToDrop.size(); ++i ) cout << _pdgToDrop[i] << ",";
      cout << endl;
    }

    if( _pdgToKeep.size()>0 ) {
      cout << "Keep these particles in the StackingAction: ";
      for( size_t i=0; i<_pdgToKeep.size(); ++i ) cout << _pdgToKeep[i] << ",";
      cout << endl;
    }

  }

  StackingAction::~StackingAction(){
  }

  void StackingAction::beginRun( double dirtG4Ymin, double dirtG4Ymax ){

    // Y limits of the dirt volume
    _dirtG4Ymin = dirtG4Ymin;
    _dirtG4Ymax = dirtG4Ymax;

    // Find the addresses of some physical volumes of interest.  See note 3.
    _dirtBodyPhysVol = G4PhysicalVolumeStore::GetInstance ()->GetVolume("DirtBody");
    _dirtCapPhysVol  = G4PhysicalVolumeStore::GetInstance ()->GetVolume("DirtCap");

  }

  // A utility function used by ClassifyNewTrack
  bool idInList( G4Track const * track, std::vector<int> const& v){
    int id(track->GetDefinition()->GetPDGEncoding());
    for( size_t i=0; i<v.size(); ++i ) {
      if( v[i] == id ) {
        return true;
      }
    }
    return false;
  }

  G4ClassificationOfNewTrack
  StackingAction::ClassifyNewTrack(const G4Track* trk){
    ++_ncalls;

    // When we get a few different tests, we can decide how to
    // make a more elegant way to call them all.
    if ( _doCosmicKiller ){
      if ( cosmicKiller(trk) ){
        return fKill;
      }
    }

    if ( !_pdgToDrop.empty() ){
      if ( idInList(trk,_pdgToDrop) ){
        return fKill;
      }
    }

    if ( !_pdgToKeep.empty() ){
      if ( !idInList(trk,_pdgToKeep) ){
        return fKill;
      }
    }

    if ( _primaryOnly ){
      if ( trk->GetParentID() != 0 ) {
        return fKill;
      }
    }

    // See Note 1) in the header file.
    if ( _killLowKineticEnergy ){
      if ( trk->GetKineticEnergy() < _eKineMin ) {
        return fKill;
      }
    }

    return fUrgent;
  }

  void StackingAction::NewStage(){
  }

  void StackingAction::PrepareNewEvent(){
    _ncalls = 0;
    ++_nevents;
  }

  bool StackingAction::cosmicKiller( const G4Track* trk){

    // Get some properties of the tracks.
    G4VPhysicalVolume* pvol = trk->GetVolume();
    G4ParticleDefinition* pdef = trk->GetDefinition();

    bool killit = false;

    if ( _killLevel > 1) {
      int pType = pdef->GetPDGEncoding();

      // just kill anything electromagnetic in the dirt
      killit =  ( pvol == _dirtBodyPhysVol &&
                  (pType == PDGCode::e_minus ||
                   pType ==  PDGCode::e_plus ||
                   PDGCode::gamma ) );
    } else {

      const G4ThreeVector& ppos = trk->GetPosition();
      G4ThreeVector p3mom = trk->GetMomentum();

      // Magic numbers for illustrative purposes.
      // You can probably get much more aggressive than this.
      // Get ycut from geometry and pcut from run time config.

      // Decide if we want to kill the track.
      killit = ( ppos.y() > (_dirtG4Ymin + _yaboveDirtYmin) && p3mom.mag() <  _cosmicpcut );
    }


    // Printout about the decision.
    if ( _nevents < 20 ) {

      // See note 4.
      G4String volName = (pvol !=0) ? pvol->GetName(): "Unknown Volume";
      G4String partName = (pdef !=0) ?
        pdef->GetParticleName() : "Unknown Particle";
      G4String killString = killit ? "Kill it": "";
      const G4ThreeVector& ppos = trk->GetPosition();
      G4ThreeVector p3mom = trk->GetMomentum();

      cout << "Cosmic Killer: "
           << setw(4)  << _nevents << " "
           << setw(4)  << _ncalls << " "
           << setw(4)  << trk->GetTrackID() << " "
           << setw(4)  << trk->GetParentID() <<  " "
           << setw(8)  << partName << "   |   "
           << ppos     << " "
           << volName  << " "
           << p3mom.mag() << "      "
           << killString
           << endl;

      // Check to see if we are in the dirt.
      if ( pvol == _dirtBodyPhysVol ){
        string tag = ( ppos.y() >= _dirtG4Ymin && ppos.y() <= _dirtG4Ymax ) ?
          "OK" : "Fail y Check";
        cout << "Cosmic Killer: tag Dirt Body " << tag << endl;
      } else if ( pvol == _dirtCapPhysVol ) {
        cout << "Cosmic Killer: tag Dirt Cap" << endl;
      }

    }

    return killit;

  } // end StackingAction::cosmicKiller

} // end namespace mu2e

