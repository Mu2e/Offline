//
// Steering routine for user tracking actions.
// If Mu2e needs many different user tracking actions, they
// should be called from this class.
//
// $Id: TrackingAction.cc,v 1.42 2014/01/21 06:22:33 kutschke Exp $
// $Author: kutschke $
// $Date: 2014/01/21 06:22:33 $
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
#include <cassert>

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Mu2eG4/inc/SteppingAction.hh"
#include "Mu2eG4/inc/TrackingAction.hh"
#include "Mu2eG4/inc/UserTrackInformation.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/ProcessCode.hh"
#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"

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
    _trajectories(nullptr),
    _sizeLimit(config.getInt("g4.particlesSizeLimit",0)),
    _currentSize(0),
    _overflowSimParticles(false),
    _mcTrajectoryMomentumCut(config.getDouble("g4.mcTrajectoryMomentumCut", 50.)),
    _saveTrajectoryMomentumCut(config.getDouble("g4.saveTrajectoryMomentumCut", 50.)),
    _mcTrajectoryMinSteps(config.getInt("g4.mcTrajectoryMinSteps", 5)),
    _steppingAction(steppingAction),
    _processInfo(0),
    _spHelper(),
    _primaryHelper()
  {

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

    // Create a user track information object and attach it to the track.
    G4VUserTrackInformation* tui = trk->GetUserInformation();
    UserTrackInformation* ti =  new UserTrackInformation();
    if (tui) {
      // G4cout << __func__ 
      // 	     << " the track was labeled as " << tui->GetType() << G4endl;
      ProcessCode cCode = Mu2eG4UserHelpers::findCreationCode(trk);
      if (cCode == ProcessCode(ProcessCode::muMinusCaptureAtRest)) {
	ti->setMuCapCode(ProcessCode::findByName((tui->GetType()).c_str()));
	// G4cout << __func__ << " set UserTrackInformation  muCapCode " 
	//        << ti->muCapCode()  << G4endl;
      }
    } 
    // else {
    //   G4cout << __func__ 
    // 	     << " the track was not labeled" << G4endl;
    // }

    // Need to cast away const-ness to do this.
    const_cast<G4Track*>(trk)->SetUserInformation(ti);

    // saveSimParticle must be called before controlTrajectorySaving.
    // but after attaching the  user track information
    saveSimParticleStart(trk);
    Mu2eG4UserHelpers::controlTrajectorySaving(trk, _sizeLimit, _currentSize, _saveTrajectoryMomentumCut);

    _steppingAction->BeginOfTrack();

    if ( !_debugList.inList() ) return;
    Mu2eG4UserHelpers::printTrackInfo( trk, "Start new Track: ", _transientMap, _timer, _mu2eOrigin);

    _timer.reset();
    _timer.start();

  }

  void TrackingAction::PostUserTrackingAction(const G4Track* trk){

    // This is safe even if it was never started.
    _timer.stop();

    // Finalize the SimParticle
    saveSimParticleEnd(trk);

    // If this particle passes the cuts, add the trajectory information to the data product.
    swapTrajectory(trk);

    // Any other clean up.
    _steppingAction->EndOfTrack();

    if ( !_debugList.inList() ) return;
    Mu2eG4UserHelpers::printTrackInfo( trk, "End Track:       ",  _transientMap, _timer, _mu2eOrigin, true);

  }

  namespace { // to use compressSimParticleCollection
    struct KeepAll {
      bool operator[](cet::map_vector_key ) const { return true; }
    };
  }
  void TrackingAction::beginEvent( const art::Handle<SimParticleCollection>& inputSimHandle,
                                   const SimParticleHelper& spHelper,
                                   const SimParticlePrimaryHelper& primaryHelper,
                                   MCTrajectoryCollection&  trajectories
                                   ) {
    _currentSize          = 0;
    _overflowSimParticles = false;
    _spHelper             = &spHelper;
    _primaryHelper        = &primaryHelper;
    _trajectories         = &trajectories;

    if(inputSimHandle.isValid()) {
      // We do not compress anything here, but use the call to reseat the pointers
      // while copying the inputs to _transientMap.
      compressSimParticleCollection(_spHelper->productID(),
                                    _spHelper->productGetter(),
                                    *inputSimHandle,
                                    KeepAll(),
                                    _transientMap);
    }
  }

  void TrackingAction::endEvent(SimParticleCollection& persistentSims ){
    Mu2eG4UserHelpers::checkCrossReferences(true,true,_transientMap);
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

    const key_type kid = _spHelper->particleKeyFromG4TrackID(trk->GetTrackID());
    const int parentId = trk->GetParentID();

    art::Ptr<GenParticle> genPtr;
    art::Ptr<SimParticle> parentPtr;

    if(parentId == 0) { // primary
      genPtr = _primaryHelper->genParticlePtr(trk->GetTrackID());
      parentPtr = _primaryHelper->simParticlePrimaryPtr(trk->GetTrackID());
    }
    else { // not a primary
      parentPtr = _spHelper->particlePtrFromG4TrackID(parentId);
    }

    // Find the physics process that created this track.
    ProcessCode creationCode = Mu2eG4UserHelpers::findCreationCode(trk);
    // we shall replace creationCode with muCapCode from UserTrackInformation if needed/present
    if (creationCode==ProcessCode(ProcessCode::muMinusCaptureAtRest)) {

      // G4cout << __func__ 
      // 	     << " particle created by " << creationCode.name()
      // 	     << " will try to replace the creation code "
      // 	     << G4endl;

      // G4VUserTrackInformation* tui = trk->GetUserInformation();
      // if (tui) {
      // 	G4cout << __func__ 
      // 	       << " the track is labeled as " << tui->GetType() 
      // 	       << G4endl;
      // 	G4cout << __func__ 
      // 	       << " muCapCode is: " 
      // 	       << (static_cast<UserTrackInformation*>(tui))->muCapCode()
      // 	       << G4endl;
      // }

      ProcessCode utic = 
	(static_cast<UserTrackInformation*>(trk->GetUserInformation()))->muCapCode();
      if (utic!=ProcessCode(ProcessCode::unknown)) {
	creationCode=utic;
      }
    }
    // G4cout << __func__ 
    // 	   << " saving particle as created by " << creationCode.name()
    // 	   << G4endl;

    // Track should not yet be in the map.  Add a debug clause to skip this test?
    if ( _transientMap.find(kid) != _transientMap.end() ){
      throw cet::exception("RANGE")
        << "SimParticle already in the event.  This should never happen. id is: "
        << kid
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
    if ( parentPtr.isNonnull() ){
      map_type::iterator i(_transientMap.find(SimParticleCollection::key_type(parentPtr.key())));
      if ( i == _transientMap.end() ){
        throw cet::exception("RANGE")
          << "Could not find parent SimParticle in PreUserTrackingAction.  id: "
          << parentPtr.key()
          << "\n";
      }
      i->second.addDaughter(_spHelper->particlePtr(trk));
    }
  }

  // Append end of track information to the existing SimParticle.
  void TrackingAction::saveSimParticleEnd(const G4Track* trk){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) return;

    key_type kid(_spHelper->particleKeyFromG4TrackID(trk->GetTrackID()));

    // Find the particle in the map.
    map_type::iterator i(_transientMap.find(kid));
    if ( i == _transientMap.end() ){
      throw cet::exception("RANGE")
        << "Could not find existing SimParticle in PostUserTrackingAction::saveSimParticleEnd()  id: "
        << kid
        << "\n";
    }

    // Reason why tracking stopped, decay, range out, etc.
    G4String pname  = Mu2eG4UserHelpers::findTrackStoppingProcessName(trk);
    ProcessCode stoppingCode(_processInfo->findAndCount(pname));

    // G4cout << __func__ 
    // 	   << " stopping process pname is " << pname << G4endl;

    // if ( pname == "muMinusCaptureAtRest") {

    //   G4VUserTrackInformation* tui = trk->GetUserInformation();
    //   if (tui) {
    // 	G4cout << __func__ 
    // 	       << " the track is labeled as " << tui->GetType() << G4endl;
    //   }

    // }

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

  // If the track passes the cuts needed to store the trajectory object, then store
  // it in the output data product.  For efficiency, the store uses a swap.
  void TrackingAction::swapTrajectory(const G4Track* trk){

    key_type kid(_spHelper->particleKeyFromG4TrackID(trk->GetTrackID()));

    std::vector<CLHEP::HepLorentzVector> const& trajectory = _steppingAction->trajectory();
    if ( int(trajectory.size()) < _mcTrajectoryMinSteps ) return;

    // Find the particle in the map.
    map_type::iterator i(_transientMap.find(kid));
    if ( i == _transientMap.end() ){
      throw cet::exception("RANGE")
        << "Could not find existing SimParticle in TrackingAction::addTrajectory  id: "
        << kid
        << "\n";
    }

    CLHEP::HepLorentzVector const& p0 = i->second.startMomentum();
    if ( p0.vect().mag() < _mcTrajectoryMomentumCut ) return;

    art::Ptr<SimParticle> sim = _spHelper->particlePtr(trk);

    // Default construct the trajectory object in the output data product.
    auto retval = _trajectories->insert( MCTrajectoryCollection::value_type( sim, MCTrajectory(sim) ));

    if ( !retval.second ){
      throw cet::exception("RANGE")
        << "In TrackingAction::addTrajectory the MCTrajectory was already present for id: "
        << kid
        << "\n";
    }

    // The data product takes ownership of the array of points that was created in SteppingAction.
    // This leaves SteppingAction with an empty array.
    MCTrajectory& traj = retval.first->second;
    _steppingAction->swapTrajectory( traj.points() );

    // So far the trajectory holds the starting point of each step.
    // Add the end point of the last step.
    traj.points().emplace_back( trk->GetPosition()-_mu2eOrigin, trk->GetGlobalTime() );

  }

} // end namespace mu2e
