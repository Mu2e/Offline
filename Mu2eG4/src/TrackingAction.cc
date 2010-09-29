//
// Steering routine for user tracking actions. 
// If Mu2e needs many different user tracking actions, they
// should be called from this class.
//
// $Id: TrackingAction.cc,v 1.7 2010/09/29 19:37:58 logash Exp $
// $Author: logash $
// $Date: 2010/09/29 19:37:58 $
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

  TrackingAction::TrackingAction( const SimpleConfig& config,
				  SteppingAction *stepping_action ):
    _debugList(),
    _physVolHelper(0),
    _currentSize(0),
    _timer(){

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

  void TrackingAction::beginEvent(){
       _spmap.clear();
       _currentSize=0;
  }


  void TrackingAction::endEvent( SimParticleCollection& simParticles ){
    saveSimParticleCopy(simParticles);
  }

  // Save start of track info to the transient store.
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
    uint32_t id       = trk->GetTrackID()-1;
    int32_t parentId  = trk->GetParentID()-1;

    // Indices into the GenParticleCollection are also 0 based.
    int32_t generatorIndex = ( parentId == -1 ) ? int32_t(id): -1;

    // Track should not yet be in the map.  Add a debug clause to skip this test?
    if ( _spmap.find(id) != _spmap.end() ){
      throw cms::Exception("RANGE")
        << "SimParticle already in the event.  This should never happen. id is: "
        << id
        << "\n";
    }

    // Add this track to the transient data.
    CLHEP::HepLorentzVector p4(trk->GetMomentum(),trk->GetTotalEnergy());
    _spmap.insert(std::make_pair(id,SimParticle( id,
                                                 parentId,
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
      std::map<uint32_t,SimParticle>::iterator i(_spmap.find(parentId));
      if ( i == _spmap.end() ){
        throw cms::Exception("RANGE")
          << "Could not find parent SimParticle in PreUserTrackingAction.  id: "
          << id
          << "\n";
      } 
      i->second.addDaughter(id);
    }

  }

  // Append end of track information to the transient store.
  void TrackingAction::saveSimParticleEnd(const G4Track* trk){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) return;

    // Persistent info uses in 0-based indices; G4 uses 1-based indices.
    uint32_t id       = trk->GetTrackID()-1;

    // Find the particle in the map.
    std::map<uint32_t,SimParticle>::iterator i(_spmap.find(id));
    if ( i == _spmap.end() ){
      throw cms::Exception("RANGE")
        << "Could not find existing SimParticle in PostUserTrackingAction.  id: "
        << id
        << "\n";
    }
    SimParticle& particle = i->second;

    // Add info about the end of the track.
    CLHEP::HepLorentzVector p4(trk->GetMomentum(),trk->GetTotalEnergy());
    particle.addEndInfo( trk->GetPosition()-_mu2eOrigin,
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
      SimParticle& particle = _spmap[id];
      cout << particle.endProperTime()*CLHEP::ns <<  " | ";
      cout << particle.startGlobalTime()*CLHEP::ns <<  " ";
      cout << particle.endGlobalTime()*CLHEP::ns <<  " | ";
      cout << _timer.cpuTime() << " " 
           << _timer.realTime() 
           << endl;
    }

    cout << endl;

  }


  // Copy transient information to the event.
  void TrackingAction::saveSimParticleCopy( SimParticleCollection& simParticles ){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) {
      edm::LogWarning("G4") << "Total of " << _currentSize 
			    << " particles were generated in the event." 
			    << endl
			    << "Only " << _sizeLimit << " are saved in output collection." 
			    << endl;
      cout << "Total of " << _currentSize 
	   << " particles were generated in the event." 
	   << endl
	   << "Only " << _sizeLimit << " are saved in output collection." 
	   << endl;
    }

    // Copy transient information to the persistent form.
    simParticles.clear();
    simParticles.reserve(_spmap.size());
    for ( std::map<uint32_t,SimParticle>::iterator i=_spmap.begin(), e=_spmap.end();
          i!=e; ++i){
      simParticles.push_back(i->second);
    }

    // Debug printout.
    /*
    PhysicalVolumeInfoCollection const& volumes = _physVolHelper->persistentInfo();
    cout << "Copy check: "
         << _spmap.size() << " "
         << simParticles.size() << " "
         << _spmap[2].startPosition() << " " 
         << simParticles[2].startPosition() << " "
         << endl;

    for ( std::map<uint32_t,SimParticle>::iterator i=_spmap.begin(), e=_spmap.end();
          i!=e; ++i){

      const vector<uint32_t>& v = i->second.daughterIds();
      size_t ndau = v.size();
      cout << "Loop:" 
           << i->first << " "
           << i->second.id() << " " 
           << i->second.parentId() << " " 
           << ndau;
      if ( ndau != 0 ){
        cout << " (";
        for ( size_t i=0; i<ndau; ++i){
          cout << " " << v[i];
        }
        cout << ")";
      }
      cout << "  |   " 
           << i->second.startVolumeIndex()  << " ";
      PhysicalVolumeInfo const& pvol = volumes.at(i->second.startVolumeIndex());
      cout << pvol.name << " "
           << pvol.copyNo;
      cout << endl;
      
    }
    */

  }

} // end namespace mu2e
