// A Producer Module that runs Geant4 and adds its output to the event.
//
// Original author Rob Kutschke
//
// Notes:



// Mu2e includes
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eHallGeom/inc/Mu2eHall.hh"
#include "Mu2eG4/inc/WorldMaker.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/Mu2eStudyWorld.hh"
#include "Mu2eG4/inc/IMu2eG4Cut.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Mu2eG4/inc/exportG4PDT.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eG4/inc/ActionInitialization.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/physicsListDecider.hh"
#include "Mu2eG4/inc/preG4InitializeTasks.hh"
#include "Mu2eG4/inc/postG4InitializeTasks.hh"
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/ExtMonFNALPixelSD.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eG4/inc/generateFieldMap.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "Mu2eG4/inc/Mu2eG4TrajectoryControl.hh"
#include "Mu2eG4/inc/Mu2eG4MultiStageParameters.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/checkConfigRelics.hh"
#include "Mu2eG4/inc/GenEventBroker.hh"
#include "Mu2eG4/inc/EventStash.hh"
#include "Mu2eG4/inc/Mu2eG4MTRunManager.hh"

// Data products that will be produced by this module.
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"

// From art and its tool chain.
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"

// Geant4 includes
#include "G4Run.hh"
#include "G4Timer.hh"
#include "G4VUserPhysicsList.hh"
#include "G4ParticleHPManager.hh"
#include "G4HadronicProcessStore.hh"
#include "G4RunManagerKernel.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4Threading.hh"

// C++ includes.
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <memory>
#include <iomanip>
#include <utility>

using namespace std;

namespace mu2e {

    class Mu2eG4 : public art::EDProducer {

public:
        Mu2eG4(fhicl::ParameterSet const& pSet);
        // Accept compiler supplied d'tor

        virtual void produce(art::Event& e) override;

        virtual void endJob() override;

        virtual void beginRun(art::Run &r) override;
        virtual void endRun(art::Run &) override;

        virtual void beginSubRun(art::SubRun &sr) override;

private:
        
        fhicl::ParameterSet pset_;
        Mu2eG4ResourceLimits mu2elimits_;
        Mu2eG4TrajectoryControl trajectoryControl_;
        Mu2eG4MultiStageParameters multiStagePars_;
      
        // The THREE functions that call new G4RunManger functions and break G4's BeamOn() into 3 pieces
        void BeamOnBeginRun( unsigned int runNumber);
        void BeamOnDoOneArtEvent( int eventNumber, G4int, const char* macroFile=0, G4int n_select=-1 );
        void BeamOnEndRun();
        
        unique_ptr<G4RunManager> _runManager;
        
        const G4int numberOfEventsToBeProcessed;
        const bool _use_G4MT;
        const G4int _nThreads;
        
        // Do we issue warnings about multiple runs?
        bool _warnEveryNewRun;

        // Do we want to export the G4 particle data table.
        bool  _exportPDTStart;
        bool  _exportPDTEnd;
      
        
        ActionInitialization const * _actionInit;
        
        //these cut objects are used in the master thread to indicate what data product is produced
        //additional thread-local cut objects are owned by ActionInitialization
        unique_ptr<IMu2eG4Cut> stackingCuts_;
        unique_ptr<IMu2eG4Cut> steppingCuts_;
        unique_ptr<IMu2eG4Cut> commonCuts_;
      
        int _rmvlevel;
        int _tmvlevel;
        int _checkFieldMap;


        // Name of a macro file to be used for controling G4 parameters after
        // the initialization phase.
        string _g4Macro;

        art::InputTag _generatorModuleLabel;

        // Helps with indexology related to persisting G4 volume information.
        // string to ptr maps, speed optimization
        // if in MT mode, only allow lookup, don't allow add
        // do a counter that counts how mnay times it was called with an unknown process
        PhysicalVolumeHelper _physVolHelper;
        
        ExtMonFNALPixelSD       *_extMonFNALPixelSD;
      
        // handles per-thread objects
        GenEventBroker _genEventBroker;
      
        // Instance name of the timeVD StepPointMC data product.
        const StepInstanceName _tvdOutputName;
        std::vector<double> timeVDtimes_;

        // Do the G4 initialization that must be done only once per job, not once per run
        void initializeG4( GeometryService& geom, art::Run const& run );

        unique_ptr<G4Timer> _timer; // local Mu2e per Geant4 event timer
        // Counters for cumulative time spent processing events by Geant4
        G4double _realElapsed;
        G4double _systemElapsed;
        G4double _userElapsed;

        const bool standardMu2eDetector_;
        G4ThreeVector originInWorld;
        
        std::vector< SensitiveDetectorHelper > SensitiveDetectorHelpers;
        EventStash _StashForEventData;
        int stashInstanceToStore;
        
        //used in testing code
        //int event_counter = 0;
        
  }; // end G4 header

    
    
Mu2eG4::Mu2eG4(fhicl::ParameterSet const& pSet):
    pset_(pSet),
    mu2elimits_(pSet.get<fhicl::ParameterSet>("ResourceLimits")),
    trajectoryControl_(pSet.get<fhicl::ParameterSet>("TrajectoryControl")),
    multiStagePars_(pSet.get<fhicl::ParameterSet>("MultiStageParameters")),
    
    numberOfEventsToBeProcessed(pSet.get<int>("numberOfEventsToProcess",0)),
    _use_G4MT(pSet.get<bool>("runinMTMode",false)),
    _nThreads(pSet.get<int>("numberOfThreads",1)),
    
    _warnEveryNewRun(pSet.get<bool>("debug.warnEveryNewRun",false)),
    _exportPDTStart(pSet.get<bool>("debug.exportPDTStart",false)),
    _exportPDTEnd(pSet.get<bool>("debug.exportPDTEnd",false)),
    
    stackingCuts_(createMu2eG4Cuts(pSet.get<fhicl::ParameterSet>("Mu2eG4StackingOnlyCut", fhicl::ParameterSet()), mu2elimits_)),
    steppingCuts_(createMu2eG4Cuts(pSet.get<fhicl::ParameterSet>("Mu2eG4SteppingOnlyCut", fhicl::ParameterSet()), mu2elimits_)),
    commonCuts_(createMu2eG4Cuts(pSet.get<fhicl::ParameterSet>("Mu2eG4CommonCut", fhicl::ParameterSet()), mu2elimits_)),
    
    // FIXME:  naming of pset parameters
    _rmvlevel(pSet.get<int>("debug.diagLevel",0)),
    _tmvlevel(pSet.get<int>("debug.trackingVerbosityLevel",0)),
    _checkFieldMap(pSet.get<int>("debug.checkFieldMap",0)),
    _g4Macro(pSet.get<std::string>("g4Macro","")),
    _generatorModuleLabel(pSet.get<std::string>("generatorModuleLabel", "")),
    _physVolHelper(),
    //_printPhysicsProcessSummary(pSet.get<bool>("debug.printPhysicsProcessSummary",false)),
    _extMonFNALPixelSD(),
    _genEventBroker(_use_G4MT),
    _tvdOutputName(StepInstanceName::timeVD),
    timeVDtimes_(pSet.get<std::vector<double> >("SDConfig.TimeVD.times")),
    _timer(std::make_unique<G4Timer>()),
    _realElapsed(0.),
    _systemElapsed(0.),
    _userElapsed(0.),
    standardMu2eDetector_((art::ServiceHandle<GeometryService>())->isStandardMu2eDetector()),
    _StashForEventData(pSet),
    stashInstanceToStore(-1)
    {
        //get the right type of RunManager depending on whether we are in MT mode or not
        if (_use_G4MT) {
            _runManager.reset(new Mu2eG4MTRunManager());
        }
        else {
            _runManager.reset(new G4RunManager());
        }

        if((_generatorModuleLabel == art::InputTag()) && multiStagePars_.genInputHits().empty()) {
            throw cet::exception("CONFIG")
            << "Error: both generatorModuleLabel and genInputHits are empty - nothing to do!\n";
    }

    //we need one SDHelper for each Worker thread, plus one extra for the Master
    //in the ActionInitialization, each worker thread is given one of these SDHs to hold its SD-related data
    //the "0th" worker thread gets the "0th" element of the vector, etc
    //we give the "_nThreads" element to the Master thread through Mu2eG4World to setup the InstanceMap in the ctor of the SDH class
    //we need only one of these SDHs to declare to art the list of products that will be produced
    SensitiveDetectorHelpers.reserve(_nThreads+1);
    
    for (int i = 0; i <= _nThreads; i++) {
        SensitiveDetectorHelpers.emplace_back(pSet.get<fhicl::ParameterSet>("SDConfig", fhicl::ParameterSet()));
    }
    
    SensitiveDetectorHelpers[0].declareProducts(this);

    produces<StatusG4>();
    produces<SimParticleCollection>();

    if(!timeVDtimes_.empty()) {
      produces<StepPointMCCollection>(_tvdOutputName.name());
    }

    if(trajectoryControl_.produce()) {
      produces<MCTrajectoryCollection>();
    }

    if(multiStagePars_.multiStage()) {
      produces<SimParticleRemapping>();
    }

    //can we simplify this and directly declare the relevent products
    //rather than contructing these unneccesary object?
    stackingCuts_->declareProducts(this);
    steppingCuts_->declareProducts(this);
    commonCuts_->declareProducts(this);
    
    produces<PhysicalVolumeInfoMultiCollection,art::InSubRun>();
     
    // The string "G4Engine" is magic; see the docs for RandomNumberGenerator.
    createEngine( art::ServiceHandle<SeedService>()->getSeed(), "G4Engine");

} // end G4:G4(fhicl::ParameterSet const& pSet);

    
  // That should really be beginJob().  G4 does not care about run
  // numbers, so we could use a hardcoded 1 for that.  The problem is
  // that Mu2e GeometryService refuses to give information outside of
  // an art run.  That makes sense for alignments and such, but
  // necessitates workarounds for G4 geometry.
void Mu2eG4::beginRun( art::Run &run){

    art::ServiceHandle<GeometryService> geom;
    SimpleConfig const& config  = geom->config();
    checkConfigRelics(config);
    
    static int ncalls(0);
    ++ncalls;

    // Do the main initialization of G4; only once per job.
    if ( ncalls == 1 ) {
      initializeG4( *geom, run );
    } else {
      if ( ncalls ==2 || _warnEveryNewRun ){
        mf::LogWarning log("G4");
        log << "G4 does not change state when we cross run boundaries - hope this is OK .... ";
        if ( ncalls == 2 && !_warnEveryNewRun ){
          log << "\nThis message will not be repeated on subsequent new runs.";
        }
      }
    }


    // Tell G4 that we are starting a new run.
    BeamOnBeginRun( run.id().run() );

    // A few more things that only need to be done only once per job,
    // not once per run, but which need to be done after the call to
    // BeamOnReadyToBeginRun.
    
    //since the cuts used by the individual threads, do we need to do this?
    if ( ncalls == 1 ) {
      stackingCuts_->finishConstruction(originInWorld);
      steppingCuts_->finishConstruction(originInWorld);
      commonCuts_->finishConstruction(originInWorld);

        //can only be run in single-threaded mode
      if( _checkFieldMap>0 && !(_use_G4MT)) generateFieldMap(originInWorld,_checkFieldMap);

      if ( _exportPDTStart ) exportG4PDT( "Start:" );//once per job
    }
}

    
void Mu2eG4::initializeG4( GeometryService& geom, art::Run const& run ){
    
    //if running in MT mode, set number of threads.
    //need to downcast the ptr to the RunManager, which was defined as a ptr to G4RunManager
    if (_use_G4MT) {
        dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->SetNumberOfThreads(_nThreads);
    }
    
    if (standardMu2eDetector_) {
      geom.addWorldG4(*GeomHandle<Mu2eHall>());
    }

    if ( _rmvlevel > 0 ) {
      mf::LogInfo logInfo("GEOM");
      logInfo << "Initializing Geant4 for " << run.id()
              << " with verbosity " << _rmvlevel << endl;
      logInfo << " Configured simParticleNumberOffset = "<< multiStagePars_.simParticleNumberOffset() << endl;
    }

    
    // Create user actions and register them with G4.
    G4VUserDetectorConstruction* allMu2e;

    //as mentioned above, we give the last element to the Master thread to setup the InstanceMap in the ctor of the SDH class
    if (standardMu2eDetector_) {

        allMu2e =
            (new WorldMaker<Mu2eWorld>(std::make_unique<Mu2eWorld>(pset_, &(SensitiveDetectorHelpers[_nThreads])  ),
                                       std::make_unique<ConstructMaterials>(pset_)) );
    }
    else {

        allMu2e =
            (new WorldMaker<Mu2eStudyWorld>(std::make_unique<Mu2eStudyWorld>(pset_, &(SensitiveDetectorHelpers[_nThreads]) ),
                                            std::make_unique<ConstructMaterials>(pset_)) );
    }
    
    
    // in the non Mu2e detector we are working in the system with the
    // origin set to 0.,0.,0. and do not use geometry service for that
    originInWorld = (!standardMu2eDetector_) ? G4ThreeVector(0.0,0.0,0.0) : (GeomHandle<WorldG4>())->mu2eOriginInWorld();
    

    preG4InitializeTasks(pset_.get<fhicl::ParameterSet>("physics"));
 
    _runManager->SetVerboseLevel(_rmvlevel);

    _runManager->SetUserInitialization(allMu2e);

    G4VUserPhysicsList* pL = physicsListDecider(pset_);
    pL->SetVerboseLevel(_rmvlevel);
    
    G4ParticleHPManager::GetInstance()->SetVerboseLevel(_rmvlevel);
    
    G4HadronicProcessStore::Instance()->SetVerbose(_rmvlevel);

    _runManager->SetUserInitialization(pL);
    

    //this is where the UserActions are instantiated
    ActionInitialization* actioninit = new ActionInitialization(pset_, _extMonFNALPixelSD, SensitiveDetectorHelpers,
                                                                &_genEventBroker, &_physVolHelper,
                                                                _use_G4MT, _nThreads, originInWorld,
                                                                mu2elimits_);
    
    //in MT mode, this is where BuildForMaster is called for master thread
    // in sequential mode, this is where Build() is called for main thread
    _runManager->SetUserInitialization(actioninit);
    
      
    // setting tracking/stepping verbosity level; tracking manager
    // sets stepping verbosity level as well;
    G4RunManagerKernel const * rmk = G4RunManagerKernel::GetRunManagerKernel();
    G4TrackingManager* tm  = rmk->GetTrackingManager();
    tm->SetVerboseLevel(_tmvlevel);
    
    // Initialize G4 for this run.
    _runManager->Initialize();

    // At this point G4 geometry and physics processes have been initialized.
    // So it is safe to modify physics processes and to compute information
    // that is derived from the G4 geometry or physics processes.

    // Mu2e specific customizations that must be done after the call to Initialize.
    postG4InitializeTasks(pset_,pL);
    
} // end G4::initializeG4


void Mu2eG4::beginSubRun(art::SubRun& sr) {
    
    unique_ptr<PhysicalVolumeInfoMultiCollection> mvi(new PhysicalVolumeInfoMultiCollection());

    if(multiStagePars_.inputPhysVolumeMultiInfo()  != art::InputTag()) {
        // Copy over data from the previous simulation stages
        art::Handle<PhysicalVolumeInfoMultiCollection> ih;
        sr.getByLabel(multiStagePars_.inputPhysVolumeMultiInfo(), ih);
        mvi->reserve(1 + ih->size());
        mvi->insert(mvi->begin(), ih->cbegin(), ih->cend());
    }

    // Append info for the current stage
    mvi->emplace_back(std::make_pair(multiStagePars_.simParticleNumberOffset(), _physVolHelper.persistentSingleStageInfo()));

    sr.put(std::move(mvi));
}


// Create one G4 event and copy its output to the art::event.
void Mu2eG4::produce(art::Event& event) {
    
    //confirm that IF we are running in MT mode we do not have inputs from previous simulation stages
    //otherwsie, throw an exception
    if (_use_G4MT) {
        
        if (   multiStagePars_.inputSimParticles() != art::InputTag()
            || multiStagePars_.inputMCTrajectories() != art::InputTag()
            || !(multiStagePars_.genInputHits().empty()) ) {
        
            throw cet::exception("CONFIG")
            << "Error: You are trying to run in MT mode with input from previous stages.  This is an invalid configuration!\n";
        }
    }
    
    //event_counter++;
    
    art::Handle<GenParticleCollection> gensHandle;
    if(!(_generatorModuleLabel == art::InputTag())) {
        event.getByLabel(_generatorModuleLabel, gensHandle);
    }

    // ProductID and ProductGetter for the SimParticleCollection.
    art::ProductID simPartId(getProductID<SimParticleCollection>(event));
    art::EDProductGetter const* simProductGetter = event.productGetter(simPartId);
    
    //stash is empty, we need to simulate events
    if (_StashForEventData.getStashSize() == 0)
    {
        
        stashInstanceToStore = 0;

        //******** these are per art::event quantities ********
        // StepPointMCCollection of input hits from the previous simulation stage
        HitHandles genInputHits;
        for(const auto& i : multiStagePars_.genInputHits()) {
            genInputHits.emplace_back(event.getValidHandle<StepPointMCCollection>(i));
        }
   
        _genEventBroker.loadEvent(genInputHits, simPartId, &event, _generatorModuleLabel, &_StashForEventData, simProductGetter);
    

        if (_use_G4MT)//in MT mode, stash size is given by the size of input GenParticleCollection
        {
            //getStashSize() can only be called after loadEvent is called
            //_StashForEventData.initializeStash(_genEventBroker.getStashSize());
            _StashForEventData.initializeStash(numberOfEventsToBeProcessed);
        }
        else//in sequential mode, the stash size is 1
        {
            _StashForEventData.initializeStash(1);
        }
        
        // Run G4 for this event and access the completed event.
        BeamOnDoOneArtEvent( event.id().event(), numberOfEventsToBeProcessed );
        
        _genEventBroker.setEventPtrToZero();

        
    }//end if stash is empty, simulate events
    
    event.put(std::move(_StashForEventData.getG4Status(stashInstanceToStore)));
    
    
    //testing stuff ********************************
//    std::cout << "in produce, printing the Stash Sim Particle info " << std::endl;
//    _StashForEventData.printInfo(stashInstanceToStore);
 
    
    //***** BEGIN HACK to reseat the SimPart Ptr, Parent Ptr, and Daughter Ptrs to point at the right place in the current art::Event
    
    std::unique_ptr<SimParticleCollection> tempSims = std::move(_StashForEventData.getSimPartCollection(stashInstanceToStore));
    
    for ( SimParticleCollection::iterator i=tempSims->begin(); i!=tempSims->end(); ++i )
    {
        //SimParticle* sim = &i->second;
        SimParticle& sim = i->second;
        
        if ( sim.isPrimary() ){
            art::Ptr<GenParticle> reseat(gensHandle, sim.genParticle().key());
            sim.genParticle() = reseat;
        }
        
        sim.parent() = art::Ptr<SimParticle>(sim.parent().id(),
                                             sim.parent().key(),
                                             simProductGetter );
        
        //the following is copied from MixMCEvents_module.cc
        std::vector<art::Ptr<SimParticle> > const& daughters = sim.daughters();
        
        if ( !daughters.empty() ) {
            std::vector<art::Ptr<SimParticle> > newDaughters;
            newDaughters.reserve(daughters.size());
            
            for ( size_t i=0; i != daughters.size(); ++i){
                art::Ptr<SimParticle> const& dau = art::Ptr<SimParticle>(daughters[i].id(), daughters[i].key(),
                                                                         simProductGetter );
                newDaughters.push_back( dau );
            }
            
            sim.setDaughterPtrs( newDaughters );
        }
        
    }//for (SimParticleCollection::iterator...
    //***** END HACK to reseat SimPart Ptrs
    
    event.put(std::move(tempSims));
    
    if(!timeVDtimes_.empty()) {
        std::unique_ptr<StepPointMCCollection> tempTVD = std::move(_StashForEventData.getTVDHits(stashInstanceToStore));
        
        for ( StepPointMCCollection::iterator i=tempTVD->begin(); i!=tempTVD->end(); ++i ){
            StepPointMC& step = *i;
            
            if ( step.simParticle().isNonnull() ){
                step.simParticle() = art::Ptr<SimParticle>(step.simParticle().id(),
                                                           step.simParticle().key(),
                                                           simProductGetter );
            }
        }
        event.put(std::move(tempTVD),_StashForEventData.getTVDName(stashInstanceToStore));
    }// if !timeVDtimes_.empty()
    
    if(trajectoryControl_.produce()) {
        //get the MCTrajCollection from the Stash and create a new one to put stuff into
        std::unique_ptr<MCTrajectoryCollection> tempTrajs = std::move(_StashForEventData.getMCTrajCollection(stashInstanceToStore));
        std::unique_ptr<MCTrajectoryCollection> outTrajectory(new MCTrajectoryCollection());
        
        for ( MCTrajectoryCollection::iterator i=tempTrajs->begin(); i!=tempTrajs->end(); ++i ){
            art::Ptr<SimParticle> newParticle(i->second.sim().id(), i->second.sim().key(), simProductGetter );
            (*outTrajectory)[newParticle] = i->second;
            (*outTrajectory)[newParticle].sim() = newParticle;
        }

        event.put(std::move(outTrajectory));
        
    }// if trajectoryControl
    
    //DO I NEED TO RESEAT THE SimParticles in the Remap?  maybe no, becasue these are only produced in sequential mode where the stash isn't used
    if(multiStagePars_.multiStage()) {
        event.put(std::move(_StashForEventData.getSimParticleRemap(stashInstanceToStore)));
    }
    
    if(SensitiveDetectorHelpers[0].extMonPixelsEnabled()) {
        std::unique_ptr<ExtMonFNALSimHitCollection> tempExtMonHits = std::move((_StashForEventData.getExtMonFNALSimHitCollection(stashInstanceToStore)));
        
        for ( ExtMonFNALSimHitCollection::iterator i=tempExtMonHits->begin(); i!=tempExtMonHits->end(); ++i ){
            ExtMonFNALSimHit& hit = *i;
            
            if ( hit.simParticle().isNonnull() ){
                hit.simParticle() = art::Ptr<SimParticle>(hit.simParticle().id(), hit.simParticle().key(), simProductGetter );
            }
        }
        event.put(std::move(tempExtMonHits));
    }//if extMonPixelsEnabled
    
    
    _StashForEventData.putSensitiveDetectorData(stashInstanceToStore, event, simProductGetter);
    _StashForEventData.putCutsData(stashInstanceToStore, event, simProductGetter);
    
    //increment the instance of the EventStash to store
    stashInstanceToStore++;
    
    if (stashInstanceToStore == _StashForEventData.getStashSize()) {
        _StashForEventData.clearStash();
    }

    
}//end Mu2eG4::produce

    
// Tell G4 that this run is over.
void Mu2eG4::endRun(art::Run & run){
      
        BeamOnEndRun();
    
}
 
    
void Mu2eG4::endJob(){

    if ( _exportPDTEnd ) exportG4PDT( "End:" );

    // Yes, these are named endRun, but they are really endJob actions.
    _physVolHelper.endRun();

    
    //this PhysicsProcessInfo has been made a per-thread item
   // G4AutoLock PIlock(&processInfoMutex);
    //if ( _printPhysicsProcessSummary ){ _processInfo.endRun(); }
    
    //_processInfo.endRun();
    
    //PIlock.unlock();
}

  
 
    
    
/**************************************************************
                    FUNCTION DEFINITIONS
 **************************************************************/
    
// Do the "begin run" parts of BeamOn.
void Mu2eG4::BeamOnBeginRun( unsigned int runNumber){

        _runManager->SetRunIDCounter(runNumber);

        bool cond = _runManager->ConfirmBeamOnCondition();
        if(!cond){
            // throw here
            return;
        }

        _realElapsed   = 0.;
        _systemElapsed = 0.;
        _userElapsed   = 0.;
    
        //this would have been set by BeamOn
        //needed for RunInitialization(), called by Run::SetNumberofEventToBeProcessed
        _runManager->SetNumberOfEventsToBeProcessed(numberOfEventsToBeProcessed);
    
        _runManager->ConstructScoringWorlds();
        _runManager->RunInitialization();

}

    
// Do the "per event" part of DoEventLoop.
void Mu2eG4::BeamOnDoOneArtEvent( int eventNumber, G4int num_events, const char* macroFile, G4int n_select){
        
        if (_use_G4MT)//MT mode
        {
            //this is where the events are actually processed
            //num_events is # of G4 events processed per art event
            
            //NOTE: G4MTRunManager::WaitForEndEventLoopWorkers() is a protected function in G4 code, so we MUST have our own MTRunManager
            //in order to access this function
            dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->InitializeEventLoop(num_events,macroFile,n_select);
            dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->Mu2eG4WaitForEndEventLoopWorkers();
            
        }
        else//sequential mode
        {
            _runManager->InitializeEventLoop(num_events,macroFile,n_select);
            
            _timer->Start();
            _runManager->ProcessOneEvent(eventNumber);
            _timer->Stop();

            // Accumulate time spent in G4 for all events in this run.
            _realElapsed   += _timer->GetRealElapsed();
            _systemElapsed += _timer->GetSystemElapsed();
            _userElapsed   += _timer->GetUserElapsed();
      
            _runManager->TerminateOneEvent();
            
        }//end if

}

    
// Do the "end of run" parts of DoEventLoop and BeamOn.
void Mu2eG4::BeamOnEndRun(){
        
        //if in sequential mode
        if (!_use_G4MT) {_runManager->TerminateEventLoop(); }
        
        //for either MT or sequential mode
        _runManager->RunTermination();
        
        //if in sequential mode
//        if (!_use_G4MT)
//        {
        G4cout << "  Event processing inside ProcessOneEvent time summary" << G4endl;
        G4cout << "  User="  << _userElapsed
        << "s Real="  << _realElapsed
        << "s Sys="   << _systemElapsed
        << "s" << G4endl;
//        }
        
}
 
} // End of namespace mu2e

using mu2e::Mu2eG4;
DEFINE_ART_MODULE(Mu2eG4);
