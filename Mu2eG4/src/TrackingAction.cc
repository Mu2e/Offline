//
// Steering routine for user tracking actions. 
// If Mu2e needs many different user tracking actions, they
// should be called from this class.
//
// $Id: TrackingAction.cc,v 1.10 2010/11/09 20:14:58 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/11/09 20:14:58 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) I anticipate that this class might eventually do a lot of
//    different jobs.  Please keep Pre/Post UserTrackingAction, 
//    free of real work - they should just dispatch other
//    methods or classes that will themselves to the real work.
//
// 2) Same comment as 1 for the beginEvent and endEvent methods.

// C++ includes
#include <iostream>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eG4/inc/TrackingAction.hh"
#include "Mu2eG4/inc/SteppingAction.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "ToyDP/inc/SimParticleCollection.hh"

// G4 incldues
#include "globals.hh"
#include "G4RunManager.hh"

using namespace std;

namespace mu2e {

  typedef SimParticleCollection::key_type key_type;

  TrackingAction::TrackingAction( const SimpleConfig& config,
                                  SteppingAction *stepping_action ):
    _debugList(),
    _physVolHelper(0),
    _timer(),
    _currentSize(0){

    _stepping = stepping_action;

    string name("g4.trackingActionEventList");
    if ( config.hasName(name) ){
      vector<int> list;
      config.getVectorInt(name,list);
      _debugList.add(list);
    }

    _sizeLimit = config.getInt("g4.particlesSizeLimit",0);
  }
  
  TrackingAction::~TrackingAction(){
  }

  void TrackingAction::PreUserTrackingAction(const G4Track* trk){

    saveSimParticleStart(trk);
    _stepping->BeginOfTrack();
  
    if ( !_debugList.inList() ) return;
    printInfo( trk, "Start new Track: ");

    _timer.reset();
    _timer.start();

  }

  void TrackingAction::PostUserTrackingAction(const G4Track* trk){

    // This is safe even if it was never started.
    _timer.stop();

    saveSimParticleEnd(trk);
    _stepping->EndOfTrack();

    if ( !_debugList.inList() ) return;
    printInfo( trk, "End Track:       ", true);

  }

  void TrackingAction::beginEvent(SimParticleCollection& sims){
    _spmap = &sims;
    _currentSize=0;
  }

  void TrackingAction::endEvent(){
    if ( !_debugList.inList() ) return;
    checkCrossReferences(true,false);
  }


  // Save start of track info.
  void TrackingAction::saveSimParticleStart(const G4Track* trk){

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        edm::LogWarning("G4") << "Maximum number of particles reached in TrackingAction: " 
                              << _currentSize << endl;
      }
      return;
    }

    // Persistent info uses in 0-based indices; G4 uses 1-based indices.
    int id       = trk->GetTrackID()-1;
    int parentId = trk->GetParentID()-1;

    // Indices into the GenParticleCollection are also 0 based.
    int32_t generatorIndex = ( parentId == -1 ) ? int32_t(id): -1;

    // Add this track; the track should not yet be in the map.
    CLHEP::HepLorentzVector p4(trk->GetMomentum(),trk->GetTotalEnergy());
    _spmap->insertOrThrow(std::make_pair(id,SimParticle( key_type(id),
                                                         key_type(parentId),
                                                         trk->GetDefinition()->GetPDGEncoding(),
                                                         generatorIndex,
                                                         trk->GetPosition()-_mu2eOrigin,
                                                         p4,
                                                         trk->GetGlobalTime(),
                                                         trk->GetProperTime(),
                                                         _physVolHelper->index(trk),
                                                         trk->GetTrackStatus(),
                                                         trk->GetWeight()
                                                         )));
    
    // If this track has a parent, tell the parent about this track.
    if ( parentId != -1 ){
      _spmap->findOrThrow(key_type(parentId)).addDaughter(key_type(id));
    }
  }

  // Append end of track information to the existing SimParticle.
  void TrackingAction::saveSimParticleEnd(const G4Track* trk){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) return;

    // Persistent info uses in 0-based indices; G4 uses 1-based indices.
    key_type id(trk->GetTrackID()-1);

    // Add info about the end of the track.  Throw if SimParticle not already there.
    CLHEP::HepLorentzVector p4(trk->GetMomentum(),trk->GetTotalEnergy());
    _spmap->findOrThrow(id).addEndInfo( trk->GetPosition()-_mu2eOrigin,
                                        p4,
                                        trk->GetGlobalTime(),
                                        trk->GetProperTime(),
                                        _physVolHelper->index(trk),
                                        trk->GetTrackStatus()
                                        );
  }

  void TrackingAction::printInfo(const G4Track* trk, const string& text, bool isEnd ){

    const G4Event* event = G4RunManager::GetRunManager()->GetCurrentEvent();

    // Get some properties of the tracks.
    G4VPhysicalVolume* pvol = trk->GetVolume();
    G4String volName = (pvol !=0) ?
      pvol->GetName(): "Unknown Volume";

    G4ParticleDefinition* pdef = trk->GetDefinition();
    G4String partName = (pdef !=0) ?
      pdef->GetParticleName() : "Unknown Particle";

    // Persistent info uses in 0-based indices; G4 uses 1-based indices.
    // Use the 0-based indices in the printout.
    uint32_t id       = trk->GetTrackID()-1;
    int32_t parentId  = trk->GetParentID()-1;

    cout << text
         << setw(5) << event->GetEventID()  << " "
         << setw(4) << id                   << " "
         << setw(4) << parentId             << " "
         << setw(8) << partName             << " | "
         << trk->GetPosition()-_mu2eOrigin  << " "
         << trk->GetMomentum()              << " "
         << trk->GetKineticEnergy()         << " "
         << volName                         << " ";

    if ( isEnd ){
      SimParticle const& particle = _spmap->at(key_type(id));
      cout << particle.endProperTime()*CLHEP::ns <<  " | ";
      cout << particle.startGlobalTime()*CLHEP::ns <<  " ";
      cout << particle.endGlobalTime()*CLHEP::ns <<  " | ";
      cout << _timer.cpuTime() << " " 
           << _timer.realTime() 
           << endl;
    }

    cout << endl;

  }

  bool TrackingAction::checkCrossReferences( bool doPrint, bool doThrow ){

    // Start by assuming we are ok; any error will turn this to false.
    bool ok(true);

    // Loop over all simulated particles.
    for ( SimParticleCollection::MapType::const_iterator i=_spmap->begin();
          i!=_spmap->end(); ++i ){

      // The next particle to look at.
      SimParticle const& sim = i->second;

      // Check that daughters point to the mother.
      std::vector<key_type> const& dau = sim.daughterIds();
      for ( size_t j=0; j<dau.size(); ++j ){
        key_type parentId = _spmap->at(key_type(dau[j])).parentId();
        if ( parentId != sim.id() ){
          
          // Daughter does not point back to the parent.
          ok = false;
          if ( doPrint ){
            edm::LogError("G4") 
              << "TrackingAction::checkCrossReferences: daughter does not point back to mother.\n";
          }
          if ( doThrow ){
            throw cms::Exception("MU2EG4") 
              << "TrackingAction::checkCrossReferences: daughter does not point back to mother.\n";
          }
        }
      }

      // Check that this particle is in the list of its parent's daughters.
      if ( sim.hasParent() ){
        key_type parentId = sim.parentId();

        // Find all daughters of this particle's mother.
        std::vector<key_type> const& mdau = _spmap->at(parentId).daughterIds();
        bool inList(false);
        for ( size_t j=0; j<mdau.size(); ++j ){
          if ( key_type(mdau[j]) == sim.id() ){
            inList = true;
            break;
          }
        }
        if ( !inList ){
          ok = false;
          if ( doPrint ){
            edm::LogError("G4") 
              << "TrackingAction::checkCrossReferences: daughter is not found amoung mother's daughters.\n";
          }
          if ( doThrow ){
            throw cms::Exception("MU2EG4") 
              << "TrackingAction::checkCrossReferences: daughter is not found amoung mother's daughters.\n";
          }
        }
      }
    }
    return ok;

  }

} // end namespace mu2e

