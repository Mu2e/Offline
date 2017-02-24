//
// Steering routine for user tracking actions.
// If Mu2e needs many different user tracking actions, they
// should be called from this class.
//
// $Id: StudyTrackingAction.cc,v 1.9 2014/08/29 22:37:41 genser Exp $
// $Author: genser $
// $Date: 2014/08/29 22:37:41 $
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
#include "Mu2eG4/inc/StudySteppingAction.hh"
#include "Mu2eG4/inc/StudyTrackingAction.hh"
#include "Mu2eG4/inc/UserTrackInformation.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/ProcessCode.hh"

// Framework includes
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"


// G4 includes
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"

using namespace std;

namespace mu2e {

  StudyTrackingAction::StudyTrackingAction( const SimpleConfig& config,
                                            StudySteppingAction* steppingAction ):
    _debugList(),
    _physVolHelper(0),
    _timer(),
    _sizeLimit(config.getInt("g4.particlesSizeLimit",0)),
    _currentSize(0),
    _overflowSimParticles(false),
    _saveTrajectoryMomentumCut(config.getDouble("g4.saveTrajectoryMomentumCut", 0.)),
    _steppingAction(steppingAction),
    _processInfo(0),
    _printTrackTiming(config.getBool("g4.printTrackTiming",true))
  {

    string name("g4.trackingActionEventList");
    if ( config.hasName(name) ){
      vector<int> list;
      config.getVectorInt(name,list);
      _debugList.add(list);
    }

  }

  StudyTrackingAction::~StudyTrackingAction(){
  }

  // Receive information that has a lifetime of a run.
  void StudyTrackingAction::beginRun( const PhysicalVolumeHelper& physVolHelper,
                                 PhysicsProcessInfo& processInfo,
                                 CLHEP::Hep3Vector const& mu2eOrigin ){
    _physVolHelper = &physVolHelper;
    _processInfo   = &processInfo;
    _mu2eOrigin    =  mu2eOrigin;
  }

  void StudyTrackingAction::endRun(){
  }

  void StudyTrackingAction::PreUserTrackingAction(const G4Track* trk){

    G4int trackingVerbosityLevel = fpTrackingManager->GetVerboseLevel();

    // Create a user track information object and attach it to the track.
    G4VUserTrackInformation* tui = trk->GetUserInformation();
    UserTrackInformation* ti =  new UserTrackInformation();

    if (tui) {
      if ( trackingVerbosityLevel > 0 ) {
        G4cout << __func__ 
               << " the track was labeled as " << tui->GetType() << G4endl;
      }
      ProcessCode cCode = Mu2eG4UserHelpers::findCreationCode(trk);
      if (cCode == ProcessCode(ProcessCode::muMinusCaptureAtRest)) {
	ti->setMuCapCode(ProcessCode::findByName((tui->GetType()).c_str()));
        if ( trackingVerbosityLevel > 0 ) {
          G4cout << __func__ << " set UserTrackInformation  muCapCode " 
                 << ti->muCapCode()  << G4endl;
        }
      }
    } 
    else {
      if ( trackingVerbosityLevel > 0 ) {
        G4cout << __func__ 
               << " the track was not labeled" << G4endl;
      }

#if G4VERSION>4100

      // here we call G4String& GetCreatorModelName() and extract the part after _
      const G4String& creatorModelName =  trk->GetCreatorModelName();
      size_t delPosition = creatorModelName.find_last_of("_");

      if ( trackingVerbosityLevel > 0 ) {
        G4cout << __func__ 
               << " full creatorModelName " 
               << creatorModelName << G4endl;
      }

      if (delPosition != G4String::npos ) {

        string modelName = creatorModelName.substr(delPosition+1);
        ProcessCode cCode = Mu2eG4UserHelpers::findCreationCode(trk);

        if ( trackingVerbosityLevel > 0 ) {
          G4cout << __func__ 
                 << " Mu2e used model name: " 
                 << modelName << G4endl;

          G4cout << __func__ 
                 << " Creator Process name from model: " 
                 << creatorModelName.substr(0,delPosition) << G4endl;

          G4cout << __func__ 
                 << " Creator Process name from trk:   " 
                 << ProcessCode::name(cCode) << G4endl;
        }

        // we label the track using the UserTrackInformation as above
 
        if (cCode == ProcessCode(ProcessCode::muMinusCaptureAtRest)) {
          ti->setMuCapCode(ProcessCode::findByName(modelName.c_str()));

          if ( trackingVerbosityLevel > 0 ) {
            G4cout << __func__ << " set UserTrackInformation  muCapCode " 
                   << ti->muCapCode()  << G4endl;
          }

        }

      }
#endif
    }

    // Need to cast away const-ness to do this.
    const_cast<G4Track*>(trk)->SetUserInformation(ti);

    // saveSimParticle must be called before controlTrajectorySaving.
    // but after attaching the  user track information
    saveSimParticleStart(trk);
    Mu2eG4UserHelpers::controlTrajectorySaving(trk, _sizeLimit, _currentSize, 
                                               _saveTrajectoryMomentumCut);

    _steppingAction->BeginOfTrack();

    if ( !_debugList.inList() ) return;
    Mu2eG4UserHelpers::printTrackInfo( trk, "Start new Track: ",  _transientMap,
                                       _timer, _mu2eOrigin);

    _timer.reset();
    _timer.start();

  }

  void StudyTrackingAction::PostUserTrackingAction(const G4Track* trk){

    // This is safe even if it was never started.
    _timer.stop();

    saveSimParticleEnd(trk);
    _steppingAction->EndOfTrack();

    if ( !_debugList.inList() ) return;
    Mu2eG4UserHelpers::printTrackInfo( trk, "End Track:       ", _transientMap,
                                       _timer, _mu2eOrigin, true, _printTrackTiming);

  }

  void StudyTrackingAction::beginEvent( art::Handle<GenParticleCollection> const& gensHandle, 
                                   art::ProductID const& simID, 
                                   art::Event const & event){
    _currentSize          = 0;
    _overflowSimParticles = false;
    _gensHandle           = &gensHandle;
    _simID                = simID;
    _event        = &event;

  }

  void StudyTrackingAction::endEvent(SimParticleCollection& persistentSims ){
    Mu2eG4UserHelpers::checkCrossReferences(true,true,_transientMap);
    persistentSims.insert( _transientMap.begin(), _transientMap.end() );
    _transientMap.clear();
    if ( !_debugList.inList() ) return;
  }

  // Save start of track info.
  void StudyTrackingAction::saveSimParticleStart(const G4Track* trk){

    G4int trackingVerbosityLevel = fpTrackingManager->GetVerboseLevel();

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in StudyTrackingAction: "
                              << _currentSize << endl;
        _overflowSimParticles = true;
      }
      return;
    }

    const int id       = trk->GetTrackID();
    const int parentId = trk->GetParentID();

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
    // we shall replace creationCode with muCapCode from UserTrackInformation if needed/present
    if (creationCode==ProcessCode(ProcessCode::muMinusCaptureAtRest)) {

      if ( trackingVerbosityLevel > 0 ) {
        G4cout << __func__ 
               << " particle created by " << creationCode.name()
               << " will try to replace the creation code "
               << G4endl;

        G4VUserTrackInformation* tui = trk->GetUserInformation();
        if (tui) {
          G4cout << __func__ 
                 << " the track is labeled as " << tui->GetType() 
                 << G4endl;
          G4cout << __func__ 
                 << " muCapCode is: " 
                 << (static_cast<UserTrackInformation*>(tui))->muCapCode()
                 << G4endl;
        }
      }

      ProcessCode utic = 
	(static_cast<UserTrackInformation*>(trk->GetUserInformation()))->muCapCode();
      if (utic!=ProcessCode(ProcessCode::unknown)) {
	creationCode=utic;
      }
    }
    if ( trackingVerbosityLevel > 0 ) {
      G4cout << __func__ 
             << " saving particle as created by " << creationCode.name()
             << G4endl;
    }

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
  void StudyTrackingAction::saveSimParticleEnd(const G4Track* trk){

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

} // end namespace mu2e
