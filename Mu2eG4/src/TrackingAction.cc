//
// Steering routine for user tracking actions.
// If Mu2e needs many different user tracking actions, they
// should be called from this class.
//
// $Id: TrackingAction.cc,v 1.34 2013/01/22 19:58:25 mjlee Exp $
// $Author: mjlee $
// $Date: 2013/01/22 19:58:25 $
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
//
// 3) Internally G4 numbers tracks 1...N.  An earlier version of this class
//    renumbered them 0...(N-1); this was an artifact of the
//    SimParticleCollection class being a std::vector, which starts at 0.
//    But now SimParticleCollection is a cet::map_vector, so it is no longer
//    necessary to do the renumbering.
//

// C++ includes
#include <iostream>

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Mu2eG4/inc/SteppingAction.hh"
#include "Mu2eG4/inc/TrackingAction.hh"
#include "Mu2eG4/inc/UserTrackInformation.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/ProcessCode.hh"

// Framework includes
#include "art/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// G4 includes
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"

using namespace std;

namespace mu2e {

  TrackingAction::TrackingAction( const SimpleConfig& config,
                                  SteppingAction     * steppingAction ):
    _debugList(),
    _physVolHelper(0),
    _timer(),
    _sizeLimit(config.getInt("g4.particlesSizeLimit",0)),
    _currentSize(0),
    _overflowSimParticles(false),
    _pointTrajectoryMomentumCut(config.getDouble("g4.pointTrajectoryMomentumCut", 50)),
    _steppingAction(steppingAction),
    _processInfo(0){

    string name("g4.trackingActionEventList");
    if ( config.hasName(name) ){
      vector<int> list;
      config.getVectorInt(name,list);
      _debugList.add(list);
    }


  }

  TrackingAction::~TrackingAction(){
  }

  // Receive information that has a lifetime of a run.
  void TrackingAction::beginRun( const PhysicalVolumeHelper& physVolHelper,
                                 PhysicsProcessInfo& processInfo,
                                 CLHEP::Hep3Vector const& mu2eOrigin ){
    _physVolHelper = &physVolHelper;
    _processInfo   = &processInfo;
    _mu2eOrigin    =  mu2eOrigin;
  }

  void TrackingAction::endRun(){
  }

  void TrackingAction::PreUserTrackingAction(const G4Track* trk){

    // saveSimParticle must be called before controlTrajectorySaving.
    saveSimParticleStart(trk);
    Mu2eG4UserHelpers::controlTrajectorySaving(trk, _sizeLimit, _currentSize, _pointTrajectoryMomentumCut);

    // Create a user trackinformation object and attach it to the track.
    // Need to cast away const-ness to do this.
    UserTrackInformation* ti =  new UserTrackInformation();
    ((G4Track *)trk)->SetUserInformation(ti);

    _steppingAction->BeginOfTrack();

    if ( !_debugList.inList() ) return;
    Mu2eG4UserHelpers::printTrackInfo( trk, "Start new Track: ", _transientMap, _timer, _mu2eOrigin);

    _timer.reset();
    _timer.start();

  }

  void TrackingAction::PostUserTrackingAction(const G4Track* trk){

    // This is safe even if it was never started.
    _timer.stop();

    saveSimParticleEnd(trk);
    _steppingAction->EndOfTrack();

    if ( !_debugList.inList() ) return;
    Mu2eG4UserHelpers::printTrackInfo( trk, "End Track:       ",  _transientMap, _timer, _mu2eOrigin, true);

  }

  void TrackingAction::beginEvent( art::Handle<GenParticleCollection> const& gensHandle, 
                                   art::ProductID const& simID, 
                                   art::Event const & event){
    _currentSize          = 0;
    _overflowSimParticles = false;
    _gensHandle           = &gensHandle;
    _simID                = simID;
    _event        = &event;

  }

  void TrackingAction::endEvent(SimParticleCollection& persistentSims ){
    Mu2eG4UserHelpers::checkCrossReferences(true,false,_transientMap);
    persistentSims.insert( _transientMap.begin(), _transientMap.end() );
    _transientMap.clear();
    if ( !_debugList.inList() ) return;
  }

  // Save start of track info.
  void TrackingAction::saveSimParticleStart(const G4Track* trk){

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in TrackingAction: "
                              << _currentSize << endl;
        _overflowSimParticles = true;
      }
      return;
    }

    int id       = trk->GetTrackID();
    int parentId = trk->GetParentID();

    // Also need the ID in this format.
    key_type kid(id);

    // Indices into the GenParticleCollection are 0 based.
    // The first n particles in the G4 track list are the same as the first n particles
    // from the generator.
    int generatorIndex = ( parentId == 0 ) ? id-1: -1;
    art::Ptr<GenParticle> genPtr;
    art::Ptr<SimParticle> parentPtr;
    if ( parentId == 0 ){
      genPtr = art::Ptr<GenParticle>(*_gensHandle,generatorIndex);
    } else{
      parentPtr = art::Ptr<SimParticle>( _simID, parentId, _event->productGetter(_simID));
    }

    // Find the physics process that created this track.
    ProcessCode creationCode = Mu2eG4UserHelpers::findCreationCode(trk);

    // Track should not yet be in the map.  Add a debug clause to skip this test?
    if ( _transientMap.find(kid) != _transientMap.end() ){
      throw cet::exception("RANGE")
        << "SimParticle already in the event.  This should never happen. id is: "
        << id
        << "\n";
    }

    // Add this track to the transient data.
    CLHEP::HepLorentzVector p4(trk->GetMomentum(),trk->GetTotalEnergy());
    _transientMap.insert(std::make_pair(kid,SimParticle( kid,
                                                         parentPtr,
                                                         static_cast<PDGCode::type>(trk->GetDefinition()->GetPDGEncoding()),
                                                         genPtr,
                                                         trk->GetPosition()-_mu2eOrigin,
                                                         p4,
                                                         trk->GetGlobalTime(),
                                                         trk->GetProperTime(),
                                                         _physVolHelper->index(trk),
                                                         trk->GetTrackStatus(),
                                                         creationCode,
                                                         trk->GetWeight()
                                                         )));

    // If this track has a parent, tell the parent about this track.
    if ( parentId != 0 ){
      map_type::iterator i(_transientMap.find(key_type(parentId)));
      if ( i == _transientMap.end() ){
        throw cet::exception("RANGE")
          << "Could not find parent SimParticle in PreUserTrackingAction.  id: "
          << id
          << "\n";
      }
      i->second.addDaughter(art::Ptr<SimParticle>( _simID, id, _event->productGetter(_simID)));
    }
  }

  // Append end of track information to the existing SimParticle.
  void TrackingAction::saveSimParticleEnd(const G4Track* trk){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) return;

    key_type id(trk->GetTrackID());

    // Find the particle in the map.
    map_type::iterator i(_transientMap.find(id));
    if ( i == _transientMap.end() ){
      throw cet::exception("RANGE")
        << "Could not find existing SimParticle in PostUserTrackingAction.  id: "
        << id
        << "\n";
    }

    // Reason why tracking stopped, decay, range out, etc.
    G4String pname  = Mu2eG4UserHelpers::findStoppingProcessName(trk);
    ProcessCode stoppingCode(_processInfo->findAndCount(pname));

    //Get kinetic energy at the begin of the last step
    double preLastStepKE = Mu2eG4UserHelpers::getPreLastStepKE(trk);


    //Get number od steps the track is made of
    int nSteps = Mu2eG4UserHelpers::getNSteps(trk);

    // Add info about the end of the track.  Throw if SimParticle not already there.
    i->second.addEndInfo( trk->GetPosition()-_mu2eOrigin,
                          CLHEP::HepLorentzVector(trk->GetMomentum(),trk->GetTotalEnergy()),
                          trk->GetGlobalTime(),
                          trk->GetProperTime(),
                          _physVolHelper->index(trk),
                          trk->GetTrackStatus(),
                          stoppingCode,
                          preLastStepKE,
                          nSteps
                          );

  }

} // end namespace mu2e
