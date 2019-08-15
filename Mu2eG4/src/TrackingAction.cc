//
// Steering routine for user tracking actions.
// If Mu2e needs many different user tracking actions, they
// should be called from this class.
//
// $Id: TrackingAction.cc,v 1.46 2014/08/29 22:37:41 genser Exp $
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
#include <cassert>

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Mu2eG4/inc/Mu2eG4SteppingAction.hh"
#include "Mu2eG4/inc/TrackingAction.hh"
#include "Mu2eG4/inc/UserTrackInformation.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"
#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "Mu2eG4/inc/Mu2eG4TrajectoryControl.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/ProcessCode.hh"
#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"

// Framework includes
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// G4 includes
#include "globals.hh"
#include "G4Event.hh"
#include "G4Ions.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4SteppingManager.hh"
#include "G4LossTableManager.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"

using namespace std;

namespace mu2e {

  TrackingAction::TrackingAction(const fhicl::ParameterSet& pset,
                                 Mu2eG4SteppingAction * steppingAction,
                                 unsigned stageOffset,
                                 const Mu2eG4TrajectoryControl& trajectoryControl,
                                 const Mu2eG4ResourceLimits& lim):
    _debugList(pset.get<std::vector<int> >("debug.trackingActionEventList", std::vector<int>())),
    _physVolHelper(0),
    _timer(),
    _trajectories(nullptr),
    _sizeLimit(lim.maxSimParticleCollectionSize()),
    _currentSize(0),
    _overflowSimParticles(false),
    _mcTrajectoryMomentumCut(trajectoryControl.mcTrajectoryMomentumCut()),
    _saveTrajectoryMomentumCut(trajectoryControl.saveTrajectoryMomentumCut()),
    _mcTrajectoryMinSteps(trajectoryControl.mcTrajectoryMinSteps()),
    _nKilledByFieldPropagator(0),
    _rangeToIgnore(pset.get<double>("physics.rangeToIgnore")),
    _steppingAction(steppingAction),
    _stageOffset(stageOffset),
    _processInfo(0),
    _printTrackTiming(pset.get<bool>("debug.printTrackTiming")),
    _spHelper(),
    _primaryHelper(),
    _stepLimitKillerVerbose(pset.get<bool>("debug.stepLimitKillerVerbose"))
  {

    if ( _stepLimitKillerVerbose && (G4Threading::G4GetThreadId() <= 0) ) {

      G4cout << __func__
             << " range threshold below which not to count killed tracks when killing slow electrons and protons: "
             << _rangeToIgnore
             << " mm"
             << G4endl;
    }

}
    
// Receive information that has a lifetime of a run.
void TrackingAction::beginRun(const PhysicalVolumeHelper* physVolHelper,
                              PhysicsProcessInfo* processInfo,
                              CLHEP::Hep3Vector const& mu2eOrigin ){
    _physVolHelper = physVolHelper;
    _processInfo = processInfo;
    _mu2eOrigin    =  mu2eOrigin;
}

  void TrackingAction::PreUserTrackingAction(const G4Track* trk){

    G4int trackingVerbosityLevel = fpTrackingManager->GetVerboseLevel();

    // Create a user track information object and attach it to the track.
    UserTrackInformation* ti =  new UserTrackInformation();

    // an old code used when mu2e was using a custom geant4 version
    G4VUserTrackInformation* tui = trk->GetUserInformation();
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
      if ( trackingVerbosityLevel > 1 ) {
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
    Mu2eG4UserHelpers::printTrackInfo( trk, "Start new Track: ", _transientMap, 
                                       _timer, _mu2eOrigin);

    _timer.reset();
    _timer.start();

}
    

void TrackingAction::PostUserTrackingAction(const G4Track* trk){

    // This is safe even if it was never started.
    _timer.stop();

    // Finalize the SimParticle
    saveSimParticleEnd(trk);

    // If this particle passes the cuts, add the trajectory
    // information to the data product.

    swapTrajectory(trk);

    // Any other clean up.
    _steppingAction->EndOfTrack();

    if ( !_debugList.inList() ) return;
    Mu2eG4UserHelpers::printTrackInfo( trk, "End Track:       ", _transientMap,
                                       _timer, _mu2eOrigin, true, _printTrackTiming);

}

    
namespace { // to use compressSimParticleCollection
    struct KeepAll {
      bool operator[](cet::map_vector_key ) const { return true; }
    };
}
    
    
void TrackingAction::beginEvent(const art::Handle<SimParticleCollection>& inputSimHandle,
                                const art::Handle<MCTrajectoryCollection>& inputTraj,
                                const SimParticleHelper& spHelper,
                                const SimParticlePrimaryHelper& primaryHelper,
                                MCTrajectoryCollection&  trajectories,
                                SimParticleRemapping& simsRemap) {
      
    _currentSize          = 0;
    _overflowSimParticles = false;
    _nKilledByFieldPropagator = 0;
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

      // old -> new particle remapping
      for(const auto& sim: *inputSimHandle) {
        art::ProductID oldID(inputSimHandle.id());
        auto key(sim.second.id().asUint());
        art::Ptr<SimParticle> oldSim(oldID, key, _spHelper->otherProductGetter(oldID));
        art::Ptr<SimParticle> newSim(_spHelper->productID(), key, _spHelper->productGetter());
        simsRemap[oldSim] = newSim;
      }
    }
      
    if(inputTraj.isValid()) {
      // Read trajectories from the previous simulation step,  reseat the pointers.
      for(const auto& i : *inputTraj) {
        const MCTrajectory& tr(i.second);
        art::Ptr<SimParticle> newSim(_spHelper->productID(), tr.sim().key(), _spHelper->productGetter());
        auto retval = _trajectories->insert(MCTrajectoryCollection::value_type(newSim, MCTrajectory(newSim)));
        if ( !retval.second ){
          throw cet::exception("RANGE")
            << "In TrackingAction::beginEvent(): error adding pre-simulated MCTrajectory for particle id "
            << newSim->id()
            << "\n";
        }
        MCTrajectory& newTraj = retval.first->second;
        newTraj.points() = tr.points();
      }
    }
    
}//beginEvent

    
void TrackingAction::endEvent(SimParticleCollection& persistentSims ){
    
    Mu2eG4UserHelpers::checkCrossReferences(true,true,_transientMap);
    persistentSims.insert( _transientMap.begin(), _transientMap.end() );
    _transientMap.clear();
      
    if ( !_debugList.inList() ) return;
}

    
// Save start of track info.
void TrackingAction::saveSimParticleStart(const G4Track* trk){
      
    G4int trackingVerbosityLevel = fpTrackingManager->GetVerboseLevel();

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
                 << (dynamic_cast<UserTrackInformation*>(tui))->muCapCode()
                 << G4endl;
        }
      }

      ProcessCode utic = 
        (dynamic_cast<UserTrackInformation*>(trk->GetUserInformation()))->muCapCode();
      if (utic!=ProcessCode(ProcessCode::unknown)) {
        creationCode=utic;
      }
    }
      
    if ( trackingVerbosityLevel > 0 ) {
      G4cout << __func__
             << " saving particle "
             << trk->GetParticleDefinition()->GetPDGEncoding()
             << ", "
             << trk->GetParticleDefinition()->GetParticleName()
             << " created by " << creationCode.name()
             << " with kinetic energy " << fixed << setw(8) << trk->GetKineticEnergy()
             << " with tot energy " << fixed << setw(8) << trk->GetTotalEnergy();
      if (trk->GetTouchable()) {
        G4cout << " in " << trk->GetTouchable()->GetVolume()->GetName()
               << " material " << trk->GetTouchable()->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName()
               << " material cuts index "
               << trk->GetTouchable()->GetVolume()->GetLogicalVolume()->GetMaterialCutsCouple()->GetIndex();
      }
      G4cout << G4endl; // step related info is not available at this stage
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

    PDGCode::type ppdgId = static_cast<PDGCode::type>(trk->GetDefinition()->GetPDGEncoding());

    // here we inspect unusual excited ions

    if (ppdgId>PDGCode::G4Threshold) {

      int excLevel = ppdgId%10;

      if ( (excLevel > 1) && (trackingVerbosityLevel > 0) ) {

        const G4ParticleDefinition* pDef = trk->GetDefinition();

        G4cout << __func__ 
               << " Warning: Unusual excited ion: " << ppdgId;
        pDef->DumpTable();
        G4cout << " Excitation energy: " 
               << dynamic_cast<const G4Ions*>(pDef)->GetExcitationEnergy() << G4endl
               << " produced by " << creationCode 
               << G4endl;
      }

    }
    
    _transientMap.insert(std::make_pair(kid,SimParticle( kid,
                                                         _stageOffset,
                                                         parentPtr,
                                                         ppdgId,
                                                         genPtr,
                                                         trk->GetPosition()-_mu2eOrigin,
                                                         p4,
                                                         trk->GetGlobalTime(),
                                                         trk->GetProperTime(),
                                                         _physVolHelper->index(trk),
                                                         trk->GetTrackStatus(),
                                                         creationCode)));

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
}//saveSimParticleStart

    
// Append end of track information to the existing SimParticle.
void TrackingAction::saveSimParticleEnd(const G4Track* trk){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) return;

    G4int trackingVerbosityLevel = fpTrackingManager->GetVerboseLevel();

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

    // adjust verbosity for isTrackKilledByFieldPropagator
    int trVerbosity = trackingVerbosityLevel;
    if ( _stepLimitKillerVerbose && trVerbosity < 1 ) trVerbosity = 1;

    if (pname == "Transportation" &&
      Mu2eG4UserHelpers::isTrackKilledByFieldPropagator(trk, trVerbosity)) {
      pname = G4String("FieldPropagator");
      if ( !(trk->GetDefinition()->GetPDGEncoding() == 11 || // electron & proton codes hardcoded for now
             trk->GetDefinition()->GetPDGEncoding() == 2212 ) ||
          G4LossTableManager::Instance()->
          GetRange(trk->GetDefinition(),
                   trk->GetKineticEnergy(),
                   trk->GetStep()->GetPreStepPoint()->GetMaterialCutsCouple())>_rangeToIgnore) {
        // count non electrons/protons or electrons/protons which have a range which is likely to make them travel
        ++_nKilledByFieldPropagator;
      }

    }

    ProcessCode stoppingCode(_processInfo->findAndCount(pname));

    if (trackingVerbosityLevel > 0 ) {
      G4int prec = G4cout.precision(15);
      const G4DynamicParticle*  pParticle = trk->GetDynamicParticle();
      double theKEnergy  = pParticle->GetKineticEnergy();
      const G4ThreeVector& theMomentumDirection = pParticle->GetMomentumDirection();
      UserTrackInformation* uti =
        (dynamic_cast<UserTrackInformation*>(trk->GetUserInformation()));
      G4cout << __func__ << " KE before int " << uti->GetKineticEnergy()
             << " Momentum direction before int " << uti->GetMomentumDirection()
             << G4endl;
      G4cout << __func__ << " KE            " << theKEnergy
             << " Momentum direction            " << theMomentumDirection
             << G4endl;
      G4cout.precision(prec);
    }


    //Get kinematics just before annihilation
    double endKE = Mu2eG4UserHelpers::getEndKE(trk);
    CLHEP::HepLorentzVector endMomentum =  Mu2eG4UserHelpers::getEndMomentum(trk);

    //Get number od steps the track is made of
    int nSteps = Mu2eG4UserHelpers::getNSteps(trk);

    // Add info about the end of the track.  Throw if SimParticle not already there.
    i->second.addEndInfo( trk->GetPosition()-_mu2eOrigin,
                          endMomentum,
                          trk->GetGlobalTime(),
                          trk->GetProperTime(),
                          _physVolHelper->index(trk),
                          trk->GetTrackStatus(),
                          stoppingCode,
                          endKE,
                          nSteps,
                          trk->GetTrackLength()
                          );

    if (trackingVerbosityLevel > 0) {
      G4int prec = G4cout.precision(15);
      G4cout << __func__
             << " particle "
             << i->second.pdgId() << ", "
             << trk->GetParticleDefinition()->GetParticleName()
             << " stopped by " << stoppingCode << ", " << pname
             << " totE deposit " << fixed << trk->GetStep()->GetTotalEnergyDeposit()
             << " NonIonE deposit " << fixed << trk->GetStep()->GetNonIonizingEnergyDeposit()
             << " in " << trk->GetVolume()->GetName()
             << " material " << trk->GetMaterial()->GetName()
             << " material cuts index " << trk->GetMaterialCutsCouple()->GetIndex()
             << G4endl;
      G4cout << __func__
             << " vertex KE " << trk->GetVertexKineticEnergy()
             << " vertex direction " << trk->GetVertexMomentumDirection()
             << G4endl;
      G4cout << __func__ << " track statuses: " << i->second.startG4Status()
             << ", " << i->second.endG4Status()
             << G4endl;
      G4cout << __func__
             << " step length " << trk->GetStepLength()
             << ", track length " << trk->GetTrackLength()
             << G4endl;
      G4cout.precision(prec);
    }
  }//saveSimParticleEnd
    
// If the track passes the cuts needed to store the trajectory object, then store
// it in the output data product.  For efficiency, the store uses a swap.
void TrackingAction::swapTrajectory(const G4Track* trk){

    key_type kid(_spHelper->particleKeyFromG4TrackID(trk->GetTrackID()));

    const auto& trajectory = _steppingAction->trajectory();
    if ( int(trajectory.size()) < _mcTrajectoryMinSteps ) return;

    // Find the particle in the map.
    map_type::iterator i(_transientMap.find(kid));
    if ( i == _transientMap.end() ){
      G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();

      mf::LogWarning("G4") << "TrackingAction::swapTrajectory: "
                           << "SimParticle is not found.\nprobably the SimParticleCollection exceeds its maximum allowed size."
                           << "Will not store MCTrajectory for: event "
                           << event->GetEventID()
                           << " Track: " << trk->GetTrackID() << "\n";
      return;
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
    traj.points().emplace_back( trk->GetPosition()-_mu2eOrigin, trk->GetGlobalTime(), trk->GetKineticEnergy() );
    
  }//swapTrajectory


} // end namespace mu2e
