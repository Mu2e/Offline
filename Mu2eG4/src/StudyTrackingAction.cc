//
// Steering routine for user tracking actions.
// If Mu2e needs many different user tracking actions, they
// should be called from this class.
//
// $Id: StudyTrackingAction.cc,v 1.1 2012/11/16 23:53:27 genser Exp $
// $Author: genser $
// $Date: 2012/11/16 23:53:27 $
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
#include "Mu2eG4/inc/StudySteppingAction.hh"
#include "Mu2eG4/inc/StudyTrackingAction.hh"
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

  StudyTrackingAction::StudyTrackingAction( const SimpleConfig& config,
                                            StudySteppingAction* steppingAction ):
    _debugList(),
    _physVolHelper(0),
    _timer(),
    _sizeLimit(config.getInt("g4.particlesSizeLimit",0)),
    _currentSize(0),
    _overflowSimParticles(false),
    _steppingAction(steppingAction),
    _processInfo(0){

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

    // saveSimParticle must be called before controlTrajectorySaving.
    saveSimParticleStart(trk);
    controlTrajectorySaving(trk);

    // Create a user trackinformation object and attach it to the track.
    // Need to cast away const-ness to do this.
    UserTrackInformation* ti =  new UserTrackInformation();
    ((G4Track *)trk)->SetUserInformation(ti);

    _steppingAction->BeginOfTrack();

    if ( !_debugList.inList() ) return;
    printInfo( trk, "Start new Track: ");

    _timer.reset();
    _timer.start();

  }

  void StudyTrackingAction::PostUserTrackingAction(const G4Track* trk){

    // This is safe even if it was never started.
    _timer.stop();

    saveSimParticleEnd(trk);
    _steppingAction->EndOfTrack();

    if ( !_debugList.inList() ) return;
    printInfo( trk, "End Track:       ", true);

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
    persistentSims.insert( _transientMap.begin(), _transientMap.end() );
    _transientMap.clear();
    checkCrossReferences(true,false);
    if ( !_debugList.inList() ) return;
  }

  // Save start of track info.
  void StudyTrackingAction::saveSimParticleStart(const G4Track* trk){

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in StudyTrackingAction: "
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
    ProcessCode creationCode = findCreationCode(trk);

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
    G4String pname  = findStoppingProcess(trk);
    ProcessCode stoppingCode(_processInfo->findAndCount(pname));

    //Get kinetic energy at the begin of the last step
    double preLastStepKE = getPreLastStepKE(trk);


    //Get number od steps the track is made of
    int nSteps = getNSteps(trk);

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

  // Enable/disable storing of trajectories based on several considerations
  void StudyTrackingAction::controlTrajectorySaving( const G4Track* trk){

    // Do not add the trajectory if the corresponding SimParticle is missing.
    if( _sizeLimit>0 && _currentSize>_sizeLimit ) return;

    bool keep = saveThisTrajectory(trk);
    G4TrackingManager* trkmgr = G4EventManager::GetEventManager()->GetTrackingManager();
    if ( keep ) {
      trkmgr->SetStoreTrajectory(true);
    } else{
      trkmgr->SetStoreTrajectory(false);
    }
  }

  // The per track decision on whether or not to store trajectories.
  bool StudyTrackingAction::saveThisTrajectory( const G4Track* trk ){

    // This is a guess at what might be useful.  Feel free to improve it.
    // We might want to change the momentum cut depending on which volume
    // the track starts in.
    CLHEP::Hep3Vector const& mom = trk->GetMomentum();
    bool keep = ( mom.mag() > 50.*CLHEP::MeV );

    return keep;
  }


  void StudyTrackingAction::printInfo(const G4Track* trk, const string& text, bool isEnd ){

    const G4Event* event = G4RunManager::GetRunManager()->GetCurrentEvent();

    // Get some properties of the tracks.
    G4VPhysicalVolume* pvol = trk->GetVolume();
    G4String volName = (pvol !=0) ?
      pvol->GetName(): "Unknown Volume";

    G4ParticleDefinition* pdef = trk->GetDefinition();
    G4String partName = (pdef !=0) ?
      pdef->GetParticleName() : "Unknown Particle";

    int id       = trk->GetTrackID();
    int parentId = trk->GetParentID();

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
      cout << trk->GetProperTime() <<  " | ";
      map_type::iterator i(_transientMap.find(key_type(id)));
      if ( i != _transientMap.end() ){
        SimParticle const& particle = i->second;
        cout << particle.startGlobalTime() <<  " ";
      } else {
        cout << -1. <<  " ";
      }

      cout << trk->GetGlobalTime() << " | ";
      cout << _timer.cpuTime() << " "
           << _timer.realTime()
           << endl;
    }

    cout << endl;

  }

  bool StudyTrackingAction::checkCrossReferences( bool doPrint, bool doThrow ){

    // Start by assuming we are ok; any error will turn this to false.
    bool ok(true);

    // Loop over all simulated particles.
    for ( map_type::const_iterator i=_transientMap.begin();
          i!=_transientMap.end(); ++i ){

      // The next particle to look at.
      SimParticle const& sim = i->second;

      // Check that daughters point to the mother.
      std::vector<key_type> const& dau = sim.daughterIds();
      for ( size_t j=0; j<dau.size(); ++j ){
        key_type parentId = _transientMap[(key_type(dau[j]))].parentId();
        if ( parentId != sim.id() ){

          // Daughter does not point back to the parent.
          ok = false;
          if ( doPrint ){
            mf::LogError("G4")
              << "StudyTrackingAction::checkCrossReferences: daughter does not point back to mother.\n";
          }
          if ( doThrow ){
            throw cet::exception("MU2EG4")
              << "StudyTrackingAction::checkCrossReferences: daughter does not point back to mother.\n";
          }
        }
      }

      // Check that this particle is in the list of its parent's daughters.
      if ( sim.hasParent() ){
        key_type parentId = sim.parentId();

        // Find all daughters of this particle's mother.
        std::vector<key_type> const& mdau = _transientMap[parentId].daughterIds();
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
            mf::LogError("G4")
              << "StudyTrackingAction::checkCrossReferences: daughter is not found amoung mother's daughters.\n";
          }
          if ( doThrow ){
            throw cet::exception("MU2EG4")
              << "StudyTrackingAction::checkCrossReferences: daughter is not found amoung mother's daughters.\n";
          }
        }
      }

    }

    return ok;

  }

  // Find the name of the process that stopped this track.
  G4String StudyTrackingAction::findStoppingProcess( G4Track const* track){

    // First check to see if Mu2e code killed this track.
    G4VUserTrackInformation* info = track->GetUserInformation();
    UserTrackInformation const* tinfo   = (UserTrackInformation*)info;

    if ( tinfo->isForced() ){
      return tinfo->code().name();
    }

    // Otherwise, G4 killed this track.
    G4VProcess const* process = track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep();

    return process->GetProcessName();
  }

  //Retrieve kinetic energy at the beginnig of the last step from UserTrackInfo
  double StudyTrackingAction::getPreLastStepKE( G4Track const* trk ) {

    G4VUserTrackInformation* info = trk->GetUserInformation();
    UserTrackInformation const* tinfo   = (UserTrackInformation*)info;
    return tinfo->preLastStepKE();

  }

  // Get the number of G4 steps the track is made of
  int StudyTrackingAction::getNSteps( G4Track const* trk ) {

    G4VUserTrackInformation* info = trk->GetUserInformation();
    UserTrackInformation const* tinfo   = (UserTrackInformation*)info;
    return tinfo->nSteps();

  }

  // Find the name of the code for the process that created this track.
  ProcessCode StudyTrackingAction::findCreationCode( G4Track const* trk){
    G4VProcess const* process = trk->GetCreatorProcess();

    // If there is no creator process, then the G4Track was created by PrimaryGenerator action.
    if ( process == 0 ){
      return ProcessCode::mu2ePrimary;
    }

    // Extract the name from the process and look up the code.
    string name = process->GetProcessName();
    return ProcessCode::findByName(name);
  }


} // end namespace mu2e
