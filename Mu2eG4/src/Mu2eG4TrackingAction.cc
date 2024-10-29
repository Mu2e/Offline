//
// Steering routine for user tracking actions.
// If Mu2e needs many different user tracking actions, they
// should be called from this class.
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
#include <algorithm>

// Mu2e includes
#include "Offline/Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Offline/Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4SteppingAction.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4TrackingAction.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4UserTrackInformation.hh"
#include "Offline/Mu2eG4/inc/SimParticleHelper.hh"
#include "Offline/Mu2eG4/inc/SimParticlePrimaryHelper.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4Config.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4TrajectoryControl.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/Mu2eUtilities/inc/compressSimParticleCollection.hh"

// Framework includes
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// G4 includes
#include "Geant4/globals.hh"
#include "Geant4/G4Event.hh"
#include "Geant4/G4Ions.hh"
#include "Geant4/G4IonTable.hh"
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4EventManager.hh"
#include "Geant4/G4SteppingManager.hh"
#include "Geant4/G4LossTableManager.hh"
#include "Geant4/G4MaterialCutsCouple.hh"
#include "Geant4/G4ParticleChange.hh"
#include "Geant4/G4DynamicParticle.hh"
#include "Geant4/G4Exp.hh"
#include "Geant4/G4Log.hh"
#include "Geant4/Randomize.hh"

using namespace std;

namespace mu2e {

  Mu2eG4TrackingAction::Mu2eG4TrackingAction(const Mu2eG4Config::Top& conf,
                                 Mu2eG4SteppingAction * steppingAction,
                                 Mu2eG4PerThreadStorage *pts):
    _debugList(conf.debug().trackingActionEventList()),
    _physVolHelper(0),
    perThreadObjects_(pts),
    _timer(),
    _sizeLimit(pts->ioconf.mu2elimits().maxSimParticleCollectionSize()),
    _currentSize(0),
    _overflowSimParticles(false),
    _mcTrajectoryMomentumCut(pts->ioconf.trajectoryControl().mcTrajectoryMomentumCut()),
    _saveTrajectoryMomentumCut(pts->ioconf.trajectoryControl().saveTrajectoryMomentumCut()),
    _mcTrajectoryMinSteps(pts->ioconf.trajectoryControl().mcTrajectoryMinSteps()),
    _nKilledByFieldPropagator(0),
    _numKilledTracks(0),
    _rangeToIgnore(conf.physics().rangeToIgnore()),
    _mu2elimits(pts->ioconf.mu2elimits()),
    _steppingAction(steppingAction),
    _processInfo(0),
    _printTrackTiming(conf.debug().printTrackTiming()),
    _stepLimitKillerVerbose(conf.debug().stepLimitKillerVerbose()),
    _muonPreAssignedDecayProperTime(-1.0),
    _muonMinPreAssignedDecayProperTime(-1.0),
    _muonMaxPreAssignedDecayProperTime(-1.0)
  {

    if ( _stepLimitKillerVerbose && (G4Threading::G4GetThreadId() <= 0) ) {

      G4cout << __func__
             << " range threshold below which not to count killed tracks when killing slow electrons and protons: "
             << _rangeToIgnore
             << " mm"
             << G4endl;
    }

    // validate muon preassigned proper time
    if (conf.physics().muonPreAssignedDecayProperTime(_muonPreAssignedDecayProperTime) &&
        conf.physics().muonMaxPreAssignedDecayProperTime(_muonMaxPreAssignedDecayProperTime)) {
      throw cet::exception("CONFIG")
        << "In Mu2eG4TrackingAction(): only one of  muonPreAssignedDecayProperTime or "
        << "muonMaxPreAssignedDecayProperTime can be set"
        << G4endl;
    }
    if (conf.physics().muonPreAssignedDecayProperTime(_muonPreAssignedDecayProperTime) &&
        conf.physics().muonMinPreAssignedDecayProperTime(_muonMinPreAssignedDecayProperTime)) {
      throw cet::exception("CONFIG")
        << "In Mu2eG4TrackingAction(): only one of  muonPreAssignedDecayProperTime or "
        << "muonMinPreAssignedDecayProperTime can be set"
        << G4endl;
    }
    if (conf.physics().muonMinPreAssignedDecayProperTime(_muonMinPreAssignedDecayProperTime) &&
        conf.physics().muonMaxPreAssignedDecayProperTime(_muonMaxPreAssignedDecayProperTime)) {
      if ( _muonMinPreAssignedDecayProperTime >= _muonMaxPreAssignedDecayProperTime ) {
        throw cet::exception("CONFIG")
          << "In Mu2eG4TrackingAction(): inconsistent "
          << "muonMinPreAssignedDecayProperTime muonMaxPreAssignedDecayProperTime: "
          << _muonMinPreAssignedDecayProperTime << " ns,"
          << _muonMaxPreAssignedDecayProperTime << " ns"
          << G4endl;
      }
    }

    if (conf.physics().muonPreAssignedDecayProperTime(_muonPreAssignedDecayProperTime)) {
      if (conf.debug().diagLevel()>0) {
        G4cout << __func__
               << " Setting muonPreAssignedDecayProperTime to "
               << _muonPreAssignedDecayProperTime << " ns" << G4endl;
      }
    }
    if (conf.physics().muonMaxPreAssignedDecayProperTime(_muonMaxPreAssignedDecayProperTime)) {
      if (conf.debug().diagLevel()>0) {
        G4cout << __func__
               << " Setting muonMaxPreAssignedDecayProperTime to "
               << _muonMaxPreAssignedDecayProperTime << " ns" << G4endl;
      }
    }
    if (conf.physics().muonMinPreAssignedDecayProperTime(_muonMinPreAssignedDecayProperTime)) {
      if (conf.debug().diagLevel()>0) {
        G4cout << __func__
               << " Setting muonMinPreAssignedDecayProperTime to "
               << _muonMinPreAssignedDecayProperTime << " ns" << G4endl;
      }
    }
  }

  // Receive information that has a lifetime of a run.
  void Mu2eG4TrackingAction::beginRun(const PhysicalVolumeHelper* physVolHelper,
                                PhysicsProcessInfo* processInfo,
                                CLHEP::Hep3Vector const& mu2eOrigin ){
    _physVolHelper = physVolHelper;
    _processInfo = processInfo;
    _mu2eOrigin    =  mu2eOrigin;
  }

  void Mu2eG4TrackingAction::PreUserTrackingAction(const G4Track* trk){

    G4int trackingVerbosityLevel = fpTrackingManager->GetVerboseLevel();

    // Create a user track information object and attach it to the track.
    Mu2eG4UserTrackInformation* ti =  new Mu2eG4UserTrackInformation();

    // an old code used when mu2e was using a custom geant4 version
    G4VUserTrackInformation* tui = trk->GetUserInformation();
    if (tui) {
      if ( trackingVerbosityLevel > 1 ) {
        G4cout << __func__
               << " the track was labeled as " << tui->GetType() << G4endl;
      }
      ProcessCode cCode = Mu2eG4UserHelpers::findCreationCode(trk);
      if (cCode == ProcessCode(ProcessCode::muMinusCaptureAtRest)) {
        ti->setMuCapCode(ProcessCode::findByName((tui->GetType()).c_str()));
        if ( trackingVerbosityLevel > 1 ) {
          G4cout << __func__ << " set Mu2eG4UserTrackInformation  muCapCode "
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

      if ( trackingVerbosityLevel > 1 ) {
        G4cout << __func__
               << " full creatorModelName "
               << creatorModelName << G4endl;
      }

      if (delPosition != G4String::npos ) {

        string modelName = creatorModelName.substr(delPosition+1);
        ProcessCode cCode = Mu2eG4UserHelpers::findCreationCode(trk);

        if ( trackingVerbosityLevel > 1 ) {
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

        // we label the track using the Mu2eG4UserTrackInformation as above

        if (cCode == ProcessCode(ProcessCode::muMinusCaptureAtRest)) {
          ti->setMuCapCode(ProcessCode::findByName(modelName.c_str()));

          if ( trackingVerbosityLevel > 1 ) {
            G4cout << __func__ << " set Mu2eG4UserTrackInformation  muCapCode "
                   << ti->muCapCode()  << G4endl;
          }

        }

      }
#endif
    }

    // Need to cast away const-ness to do this.
    const_cast<G4Track*>(trk)->SetUserInformation(ti);

    if (_muonPreAssignedDecayProperTime>=0.0) {
      const G4DynamicParticle* dynPart = trk->GetDynamicParticle();
      if ((dynPart->GetPDGcode() == PDGCode::mu_minus) or (dynPart->GetPDGcode() == PDGCode::mu_plus)) {
        // if track is a muon and preassigned proper time is set, assing it
        if (trackingVerbosityLevel>0) {
          G4cout << __func__
                 << " Setting muonPreAssignedDecayProperTime to track "
                 << trk->GetTrackID() << G4endl;
        }
        const_cast<G4DynamicParticle*>(dynPart)
          ->SetPreAssignedDecayProperTime(_muonPreAssignedDecayProperTime*CLHEP::ns);
      }
    }

    // bias muon proper time if requested
    if ((_muonMaxPreAssignedDecayProperTime>=0.0) ||
        (_muonMinPreAssignedDecayProperTime>=0.0)) {
      const G4DynamicParticle* dynPart = trk->GetDynamicParticle();
      if ((dynPart->GetPDGcode() == PDGCode::mu_minus) or (dynPart->GetPDGcode() == PDGCode::mu_plus)) {
        // if track is a muon and preassigned proper time is set, assing it
        if (trackingVerbosityLevel>0) {
          G4cout << __func__
                 << " Setting random muonPreAssignedDecayProperTime to the track "
                 << " between muonMaxPreAssignedDecayProperTime and muonMaxPreAssignedDecayProperTime to track "
                 << trk->GetTrackID() << G4endl;
        }
        // 0 is a special case
        if (_muonMaxPreAssignedDecayProperTime == 0.0) {
          const_cast<G4DynamicParticle*>(dynPart)
            ->SetPreAssignedDecayProperTime(_muonMaxPreAssignedDecayProperTime*CLHEP::ns);
        } else {
          double pdgLifeTime = dynPart->GetParticleDefinition()->GetPDGLifeTime();
          double expNormalizedMaxProperTime =
            (_muonMaxPreAssignedDecayProperTime>=0.0) ?
            G4Exp(-1.0*_muonMaxPreAssignedDecayProperTime*CLHEP::ns/pdgLifeTime) :
            0.0;
          double expNormalizedMinProperTime =
            (_muonMinPreAssignedDecayProperTime>=0.0) ?
            G4Exp(-1.0*_muonMinPreAssignedDecayProperTime*CLHEP::ns/pdgLifeTime) :
            1.0;
          double offset = std::max(expNormalizedMaxProperTime, 0.0);
          double width  = std::min(expNormalizedMinProperTime, 1.0) - offset;
          // transforming the value into the potentially narrower interval
          double theValue = offset + G4UniformRand()*width;
          const_cast<G4DynamicParticle*>(dynPart)
            ->SetPreAssignedDecayProperTime(-1.0*G4Log(theValue)*pdgLifeTime);
        }
      }
    }

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


  void Mu2eG4TrackingAction::PostUserTrackingAction(const G4Track* trk){

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


  void Mu2eG4TrackingAction::beginEvent() {
    const Mu2eG4IOConfigHelper& ioconf = perThreadObjects_->ioconf;

    // Read in data products from previous stages and reseat SimParticle pointers
    // The returned object is allowed to be in an invalid state for GenParticle driven jobs
    // and non-filtered events in subsequent stages which do not have any primary
    // StepPointMCs, for example.
    auto simsInfo = ioconf.inputs().inputSimParticles(*perThreadObjects_->artEvent);
    if( simsInfo.isValid()) {
      const SimParticleCollection& inputSims = simsInfo.sims.ref();
      // We do not compress anything here, but use the call to reseat the pointers
      // while copying the inputs to _transientMap.
      compressSimParticleCollection(perThreadObjects_->simParticleHelper->productID(),
                                    perThreadObjects_->simParticleHelper->productGetter(),
                                    inputSims,
                                    KeepAll(),
                                    _transientMap);

      // old -> new particle remapping
      for(const auto& sim: inputSims) {
        art::ProductID oldID(simsInfo.id);
        auto key(sim.second.id().asUint());
        art::Ptr<SimParticle> oldSim(oldID, key, perThreadObjects_->simParticleHelper->otherProductGetter(oldID));
        art::Ptr<SimParticle> newSim(perThreadObjects_->simParticleHelper->productID(), key, perThreadObjects_->simParticleHelper->productGetter());
        (*perThreadObjects_->simRemapping)[oldSim] = newSim;
      }

      if(art::InputTag() != perThreadObjects_->ioconf.inputs().inputMCTrajectories()) {
        auto const& inputTraj = perThreadObjects_->artEvent->getValidHandle<MCTrajectoryCollection>(ioconf.inputs().inputMCTrajectories());

        for(const auto& i : *inputTraj) {
          const MCTrajectory& tr(i.second);
          art::Ptr<SimParticle> newSim(perThreadObjects_->simParticleHelper->productID(), tr.sim().key(), perThreadObjects_->simParticleHelper->productGetter());
          auto retval = perThreadObjects_->mcTrajectories->insert(MCTrajectoryCollection::value_type(newSim, MCTrajectory(newSim)));
          if ( !retval.second ){
            throw cet::exception("RANGE")
              << "In Mu2eG4TrackingAction::beginEvent(): error adding pre-simulated MCTrajectory for particle id "
              << newSim->id()
              << "\n";
          }
          MCTrajectory& newTraj = retval.first->second;
          newTraj.points() = tr.points();
        }
      }
    }

    _currentSize          = 0;
    _overflowSimParticles = false;
    _nKilledByFieldPropagator = 0;
    _numKilledTracks = 0;

  }//beginEvent


  void Mu2eG4TrackingAction::endEvent(){

    Mu2eG4UserHelpers::checkCrossReferences(true,true,_transientMap);
    perThreadObjects_->simPartCollection->insert( _transientMap.begin(), _transientMap.end() );
    _transientMap.clear();

    if ( !_debugList.inList() ) return;
  }


  // Save start of track info.
  void Mu2eG4TrackingAction::saveSimParticleStart(const G4Track* trk){

    G4int trackingVerbosityLevel = fpTrackingManager->GetVerboseLevel();

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in Mu2eG4TrackingAction: "
                             << _currentSize << endl;
        _overflowSimParticles = true;
      }
      return;
    }

    const key_type kid = perThreadObjects_->simParticleHelper->particleKeyFromG4TrackID(trk->GetTrackID());

    const int parentId = trk->GetParentID();

    art::Ptr<GenParticle> genPtr;
    art::Ptr<SimParticle> parentPtr;
    ProcessCode creationCode{ProcessCode::unknown};

    if(parentId == 0) { // primary
      const auto g4TrkID = trk->GetTrackID();
      genPtr = perThreadObjects_->simParticlePrimaryHelper->genParticlePtr(g4TrkID);
      parentPtr = perThreadObjects_->simParticlePrimaryHelper->simParticlePrimaryPtr(g4TrkID);

      const auto orig = perThreadObjects_->simParticlePrimaryHelper->getEntry(g4TrkID);
      // A particle that is a continuation of a particle from an earlier stage
      // in Mu2eG4 multistage jobs gets mu2ePrimary creation code.
      // The StageParticle case is different, it is treated more like
      // a custom physics process with a set of dedicated creation codes.

      if(auto stpart = std::get_if<const StageParticle*>(&orig)) {
        creationCode = (*stpart)->creationCode();
      }
      else {
        creationCode = ProcessCode::mu2ePrimary;
      }

    }
    else { // not a primary
      parentPtr = perThreadObjects_->simParticleHelper->particlePtrFromG4TrackID(parentId);

      // Find the physics process that created this track.
      creationCode = Mu2eG4UserHelpers::findCreationCode(trk);
      // we shall replace creationCode with muCapCode from Mu2eG4UserTrackInformation if needed/present

      if (creationCode==ProcessCode(ProcessCode::muMinusCaptureAtRest)) {

        if ( trackingVerbosityLevel > 1 ) {
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
                   << (dynamic_cast<Mu2eG4UserTrackInformation*>(tui))->muCapCode()
                   << G4endl;
          }
        }

        ProcessCode utic =
          (dynamic_cast<Mu2eG4UserTrackInformation*>(trk->GetUserInformation()))->muCapCode();
        if (utic!=ProcessCode(ProcessCode::unknown)) {
          creationCode=utic;
        }
      }
    } // else - not a primary


    if ( trackingVerbosityLevel > 1 ) {
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

    // printing ions
    if ( trackingVerbosityLevel > 1 && ppdgId>PDGCode::G4Threshold ) {
      G4cout << __func__ << " Ion pdgid:          "
             << ppdgId
             << G4endl;
      const G4ParticleDefinition* pDef = trk->GetDefinition();
      pDef->DumpTable();
      // print data specific to ions
      // Get excitation energy of nucleus
      // Charge
      // Get Isomer level (=0 for ground state)
      const G4Ions* pG4Ion = dynamic_cast<const G4Ions*>(pDef);
      G4int flbi = pG4Ion->GetFloatLevelBaseIndex();
      G4cout << __func__ << " Ion specific data:  "
             << pG4Ion->GetExcitationEnergy()
             << ", " << pG4Ion->GetPDGCharge()
             << ", " << pG4Ion->GetIsomerLevel()
             << ", " << flbi
        //   << ", " << static_cast<G4int>(pG4Ion->GetFloatLevelBase())
             << ", " << std::string(1,G4Ions::FloatLevelBaseChar(G4Ions::FloatLevelBase(flbi)))
             << G4endl;
      Mu2eG4UserHelpers::printTrackInfo( trk, " Ion:          ", _transientMap,
                                         _timer, _mu2eOrigin);
    }

    // extracting info specific to ions
    SimParticle::IonDetail ion;
    if ( ppdgId>PDGCode::G4Threshold ) {
      const G4ParticleDefinition* pDef = trk->GetDefinition();
      ion.excitationEnergy = dynamic_cast<const G4Ions*>(pDef)->GetExcitationEnergy();
      ion.floatLevelBaseIndex = dynamic_cast<const G4Ions*>(pDef)->GetFloatLevelBaseIndex();
    }

    _transientMap.insert(std::make_pair(kid,SimParticle( kid,
                                                         perThreadObjects_->simParticleHelper->simStage(),
                                                         parentPtr,
                                                         ppdgId,
                                                         genPtr,
                                                         trk->GetPosition()-_mu2eOrigin,
                                                         p4,
                                                         trk->GetGlobalTime(),
                                                         trk->GetProperTime(),
                                                         _physVolHelper->index(trk),
                                                         trk->GetTrackStatus(),
                                                         creationCode,
                                                         ion)));

    // If this track has a parent, tell the parent about this track.
    if ( parentPtr.isNonnull() ){
      map_type::iterator i(_transientMap.find(SimParticleCollection::key_type(parentPtr.key())));
      if ( i == _transientMap.end() ){
        throw cet::exception("RANGE")
          << "Could not find parent SimParticle in " << __func__ << ".  id: "
          << parentPtr.key()
          << "\n";
      }
      i->second.addDaughter(perThreadObjects_->simParticleHelper->particlePtr(trk));

      // // print parent of an ion
      //
      // int parPDGId = i->second.pdgId();
      // if ( ppdgId >PDGCode::G4Threshold ) {
      //   G4String pName = "";
      //   if ( parPDGId >PDGCode::G4Threshold ) {
      //     G4cout << __func__ << " Parent base ion pdgid: " << parPDGId-parPDGId%10 << G4endl;
      //     G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(parPDGId-parPDGId%10);
      //     pName = (ion) ? ion->GetParticleName() : "Unknown";
      //   } else {
      //     pName = G4ParticleTable::GetParticleTable()->FindParticle(parPDGId)->GetParticleName();
      //   }
      //   G4cout << __func__ << " Ion parent with approximate name : "
      //          << i->second.id()
      //          << ", " << parPDGId
      //          << ", " << pName
      //          << ", created by " << i->second.creationCode().name()
      //          << ", stopped by " << i->second.stoppingCode().name()
      //          << G4endl;
      // }
      // // print if parent is an ion
      // if ( parPDGId)) {
      //   G4cout << __func__ << " Ion daughter pdgid: " << ppdgId << G4endl;
      //   Mu2eG4UserHelpers::printTrackInfo( trk, "ion daughter: ", _transientMap,
      //                                      _timer, _mu2eOrigin);
      // }
    }
  }//saveSimParticleStart


  // Append end of track information to the existing SimParticle.
  void Mu2eG4TrackingAction::saveSimParticleEnd(const G4Track* trk){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) return;

    G4int trackingVerbosityLevel = fpTrackingManager->GetVerboseLevel();

    key_type kid(perThreadObjects_->simParticleHelper->particleKeyFromG4TrackID(trk->GetTrackID()));

    // Find the particle in the map.
    map_type::iterator i(_transientMap.find(kid));
    if ( i == _transientMap.end() ){
      throw cet::exception("RANGE")
        << "Could not find existing SimParticle in Mu2eG4TrackingAction::saveSimParticleEnd()  id: "
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
      pname = G4String("mu2eFieldPropagator");
      if ( !(trk->GetDefinition()->GetPDGEncoding() == PDGCode::e_minus ||
             trk->GetDefinition()->GetPDGEncoding() == PDGCode::proton ) ||
           G4LossTableManager::Instance()->
           GetRange(trk->GetDefinition(),
                    trk->GetKineticEnergy(),
                    trk->GetStep()->GetPreStepPoint()->GetMaterialCutsCouple())>_rangeToIgnore) {
        // count non electrons/protons or electrons/protons which have a range which is likely to make them travel
        ++_nKilledByFieldPropagator;
      }

    }

    if ( pname == "mu2eSpecialCutsProcess"
         && trk->GetCurrentStepNumber() >= static_cast<int>(_mu2elimits.maxStepsPerTrack())) {
      if ( _stepLimitKillerVerbose ) {
        G4cout << __func__ << " WARNING: kill particle in "
               << trk->GetStep()->GetPreStepPoint()->
          GetPhysicalVolume()->GetLogicalVolume()->GetName()
               << " due to large number of steps." << G4endl;
        Mu2eG4UserHelpers::printKilledTrackInfo(trk);
      }
      ++_numKilledTracks;

      // Changing the pname to the old name used in such cases
      pname = G4String("mu2eMaxSteps");

    }

    ProcessCode stoppingCode(_processInfo->findAndCount(pname));


    //Get kinematics just before the post step doit process acted (e.g., decay, annihilation, etc.)
    //note that while we have the intermediate position we use it only for diagnostic printouts
    double endKE = Mu2eG4UserHelpers::getEndKE(trk);
    CLHEP::HepLorentzVector endMomentum =  Mu2eG4UserHelpers::getEndMomentum(trk);
    double endGlobalTime = Mu2eG4UserHelpers::getEndGlobalTime(trk);
    double endProperTime = Mu2eG4UserHelpers::getEndProperTime(trk);

    //Get number of steps the track is made of
    int nSteps = Mu2eG4UserHelpers::getNSteps(trk);

    if (trackingVerbosityLevel > 1 ) {
      G4int prec = G4cout.precision(15);
      const G4DynamicParticle*  pParticle = trk->GetDynamicParticle();
      double theKEnergy  = pParticle->GetKineticEnergy();
      const G4ThreeVector& theMomentumDirection = pParticle->GetMomentumDirection();
      Mu2eG4UserTrackInformation* uti =
        (dynamic_cast<Mu2eG4UserTrackInformation*>(trk->GetUserInformation()));
      G4StepPoint const* lastPreStepPoint = trk->GetStep()->GetPreStepPoint();
      G4cout << __func__ << " KE pre step   "  << lastPreStepPoint->GetKineticEnergy()
             << " Momentum direction pre step   " << lastPreStepPoint->GetMomentumDirection()
             << " Position pre step   " << lastPreStepPoint->GetPosition()
             << " Global time pre step   " << lastPreStepPoint->GetGlobalTime()
             << " Proper time pre step   " << lastPreStepPoint->GetProperTime()
             << G4endl;
      G4cout << __func__ << " KE before int " << uti->GetKineticEnergy()
             << " Momentum direction before int " << uti->GetMomentumDirection()
             << " Position before int " << uti->GetPosition()
             << " Global time before int " << uti->GetGlobalTime()
             << " Proper time before int " << uti->GetProperTime()
             << G4endl;
      G4cout << __func__ << " KE            " << theKEnergy
             << " Momentum direction            " << theMomentumDirection
             << " Position            " << trk->GetPosition()
             << " Global time            " << trk->GetGlobalTime()
             << " Proper time            " << trk->GetProperTime()
             << G4endl;
      G4cout << __func__ << " KE stored     " << endKE
             << " G4 Position stored  " << trk->GetPosition()
             << " Global time stored     " << endGlobalTime
             << " Proper time stored     " << endProperTime
             << G4endl;
      G4cout.precision(prec);
    }

    // Add info about the end of the track.  Throw if SimParticle not already there.
    i->second.addEndInfo( trk->GetPosition()-_mu2eOrigin,
                          endMomentum, // based on pre last step
                          endGlobalTime, // based on pre last step
                          endProperTime, // based on pre last step
                          _physVolHelper->index(trk),
                          trk->GetTrackStatus(),
                          stoppingCode,
                          endKE, // based on pre last step
                          nSteps,
                          trk->GetTrackLength()
                          );

    // the following is to establish for debugging if the parent pdgid to see if it is an ion
    // const int parentId = trk->GetParentID();
    // art::Ptr<SimParticle> parentPtr;
    // int parPDGId = 0;
    // if(parentId == 0) { // primary
    //   parentPtr = perThreadObjects_->simParticlePrimaryHelper->simParticlePrimaryPtr(trk->GetTrackID());
    // }
    // else { // not a primary
    //   parentPtr = perThreadObjects_->simParticleHelper->particlePtrFromG4TrackID(parentId);
    // }
    // if ( parentPtr.isNonnull() ){
    //   map_type::iterator i(_transientMap.find(SimParticleCollection::key_type(parentPtr.key())));
    //   if ( i == _transientMap.end() ){
    //     throw cet::exception("RANGE")
    //       << "Could not find parent SimParticle in " << __func__ << ".  id: "
    //       << parentPtr.key()
    //       << "\n";
    //   }
    //   parPDGId = i->second.pdgId();
    // }

    if ( trackingVerbosityLevel > 1
         // || trk->GetDefinition()->GetPDGEncoding()>PDGCode::G4Threshold
        ) {
      G4int prec = G4cout.precision(15);
      G4cout << __func__
             << " particle "
             << i->second.pdgId() << ", "
             << trk->GetParticleDefinition()->GetParticleName()
             << " stopped by " << stoppingCode // << ", " << pname
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
             << ", step number " << trk->GetCurrentStepNumber()
             << ", created by " << Mu2eG4UserHelpers::findCreationCode(trk)
             << G4endl;
      if ( trk->GetDefinition()->GetPDGEncoding()>PDGCode::G4Threshold ) {
        const G4Ions* pG4Ion = dynamic_cast<const G4Ions*>(trk->GetDefinition());
        G4int flbi = pG4Ion->GetFloatLevelBaseIndex();
        G4cout << __func__ << " Ion "
               << pG4Ion->GetParticleName()
               << " with excitaion energy: "
               << pG4Ion->GetExcitationEnergy()
               << " with float level base: "
               << std::string(1,G4Ions::FloatLevelBaseChar(G4Ions::FloatLevelBase(flbi))) << " " << flbi
               << G4endl;
      }
      G4cout.precision(prec);
    }
  }//saveSimParticleEnd

  // If the track passes the cuts needed to store the trajectory object, then store
  // it in the output data product.  For efficiency, the store uses a swap.
  void Mu2eG4TrackingAction::swapTrajectory(const G4Track* trk){

    key_type kid(perThreadObjects_->simParticleHelper->particleKeyFromG4TrackID(trk->GetTrackID()));

    const auto& trajectory = _steppingAction->trajectory();
    if ( int(trajectory.size()) < _mcTrajectoryMinSteps ) return;

    // Find the particle in the map.
    map_type::iterator i(_transientMap.find(kid));
    if ( i == _transientMap.end() ){
      G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();

      mf::LogWarning("G4") << "Mu2eG4TrackingAction::swapTrajectory: "
                           << "SimParticle is not found.\nprobably the SimParticleCollection exceeds its maximum allowed size."
                           << "Will not store MCTrajectory for: event "
                           << event->GetEventID()
                           << " Track: " << trk->GetTrackID() << "\n";
      return;
    }

    CLHEP::HepLorentzVector const& p0 = i->second.startMomentum();
    if ( p0.vect().mag() < _mcTrajectoryMomentumCut ) return;

    art::Ptr<SimParticle> sim = perThreadObjects_->simParticleHelper->particlePtr(trk);

    // Default construct the trajectory object in the output data product.
    auto retval = perThreadObjects_->mcTrajectories->insert( MCTrajectoryCollection::value_type( sim, MCTrajectory(sim) ));

    if ( !retval.second ){
      throw cet::exception("RANGE")
        << "In Mu2eG4TrackingAction::addTrajectory the MCTrajectory was already present for id: "
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
