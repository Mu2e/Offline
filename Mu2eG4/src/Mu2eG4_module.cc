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
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eG4/inc/generateFieldMap.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "Mu2eG4/inc/Mu2eG4TrajectoryControl.hh"
#include "Mu2eG4/inc/Mu2eG4MultiStageParameters.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/checkConfigRelics.hh"
#include "Mu2eG4/inc/GenEventBroker.hh"
#include "Mu2eG4/inc/Mu2eG4MTRunManager.hh"
#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined G4VIS_USE_OPENGLQT )
#include "Mu2eG4/inc/Mu2eVisCommands.hh"
#endif

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
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"

// Geant4 includes
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined G4VIS_USE_OPENGLQT )
#include "G4VisExecutive.hh"
#endif
#include "G4Run.hh"
#include "G4Timer.hh"
#include "G4VUserPhysicsList.hh"
#include "G4ParticleHPManager.hh"
#include "G4HadronicProcessStore.hh"
#include "G4RunManagerKernel.hh"
#include "G4SDManager.hh"

// C++ includes.
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <memory>
#include <iomanip>
#include <utility>

using namespace std;

namespace {
  art::InputTag const invalid_tag{};
}

namespace mu2e {

  class Mu2eG4 : public art::SharedProducer {
  public:
    Mu2eG4(fhicl::ParameterSet const& pSet, art::ProcessingFrame const& pf);
    // Accept compiler supplied d'tor

  private:
      void produce(art::Event& e, art::ProcessingFrame const& pf) override;
      void endJob(art::ProcessingFrame const& pf) override;
      void beginRun(art::Run &r, art::ProcessingFrame const& pf) override;
      void endRun(art::Run &r, art::ProcessingFrame const& pf) override;
      void beginSubRun(art::SubRun &sr, art::ProcessingFrame const& pf) override;

    fhicl::ParameterSet pset_;
    Mu2eG4ResourceLimits mu2elimits_;
    Mu2eG4MultiStageParameters multiStagePars_;

    // The THREE functions that call new G4RunManger functions and break G4's BeamOn() into 3 pieces
    void BeamOnBeginRun( unsigned int runNumber);
    void BeamOnDoOneArtEvent( int eventNumber, G4int, const char* macroFile=0, G4int n_select=-1 );
    void BeamOnEndRun();


    unique_ptr<G4RunManager> _runManager;

    const bool _use_G4MT;
    const G4int _nThreads;

    // Do we issue warnings about multiple runs?
    bool _warnEveryNewRun;

    // Do we want to export the G4 particle data table.
    bool  _exportPDTStart;
    bool  _exportPDTEnd;

    // to be able to make StorePhysicsTable call after the event loop started
    G4VUserPhysicsList* physicsList_;

    ActionInitialization const * _actionInit;

    //these cut objects are used in the master thread to indicate what data product is produced
    //additional thread-local cut objects are owned by ActionInitialization
//TAKING THESE OUT TEMPORARILY
//    unique_ptr<IMu2eG4Cut> stackingCuts_;
//    unique_ptr<IMu2eG4Cut> steppingCuts_;
//    unique_ptr<IMu2eG4Cut> commonCuts_;


    int _rmvlevel;
    int _tmvlevel;
    int _smvlevel;
    int _checkFieldMap;


<<<<<<< HEAD
        // handles per-thread objects
        GenEventBroker _genEventBroker;
=======
    art::InputTag _generatorModuleLabel;

    // Helps with indexology related to persisting G4 volume information.
    // string to ptr maps, speed optimization
    // if in MT mode, only allow lookup, don't allow add
    // do a counter that counts how mnay times it was called with an unknown process
    PhysicalVolumeHelper _physVolHelper;
>>>>>>> First stage of changes to code to integrate G4MT into art3 - get rid of stashes

    vector< SensitiveDetectorHelper > SensitiveDetectorHelpers;
    ExtMonFNALPixelSD       *_extMonFNALPixelSD;

      
    // handles per-thread objects
    GenEventBroker _genEventBroker;
    

    // Do the G4 initialization that must be done only once per job, not once per run
    void initializeG4( GeometryService& geom, art::Run const& run );

<<<<<<< HEAD
        const bool    _standardMu2eDetector;
        G4ThreeVector _originInWorld;
=======
    unique_ptr<G4Timer> _timer; // local Mu2e per Geant4 event timer
    // Counters for cumulative time spent processing events by Geant4
    G4double _realElapsed;
    G4double _systemElapsed;
    G4double _userElapsed;
>>>>>>> First stage of changes to code to integrate G4MT into art3 - get rid of stashes

    const bool standardMu2eDetector_;
    G4ThreeVector originInWorld;

    // In sequential mode, there is one SensitiveDetetorHelper object and it is
    // at index 0; in MT mode the SensitiveDetectorHelper object for the master
    // thread is at index _nThreads.
    size_t _masterThreadIndex;
      
    // count the number of events that have been excluded because they did not
    // pass the filtering in Mu2eG4EventAction
    int numExcludedEvents = 0;
      
      
    CLHEP::HepJamesRandom _engine;

  }; // end G4 header

  Mu2eG4::Mu2eG4(fhicl::ParameterSet const& pSet, art::ProcessingFrame const& procFrame):
    SharedProducer{pSet},
    pset_(pSet),
    mu2elimits_(pSet.get<fhicl::ParameterSet>("ResourceLimits")),
    multiStagePars_(pSet.get<fhicl::ParameterSet>("MultiStageParameters")),

    _use_G4MT(pSet.get<bool>("runinMTMode",false)),
    _nThreads(pSet.get<int>("numberOfThreads",1)),

    _warnEveryNewRun(pSet.get<bool>("debug.warnEveryNewRun",false)),
    _exportPDTStart(pSet.get<bool>("debug.exportPDTStart",false)),
    _exportPDTEnd(pSet.get<bool>("debug.exportPDTEnd",false)),


//    stackingCuts_(createMu2eG4Cuts(pSet.get<fhicl::ParameterSet>("Mu2eG4StackingOnlyCut", {}), mu2elimits_)),
//    steppingCuts_(createMu2eG4Cuts(pSet.get<fhicl::ParameterSet>("Mu2eG4SteppingOnlyCut", {}), mu2elimits_)),
//    commonCuts_(createMu2eG4Cuts(pSet.get<fhicl::ParameterSet>("Mu2eG4CommonCut", {}), mu2elimits_)),

    _rmvlevel(pSet.get<int>("debug.diagLevel",0)),
    _tmvlevel(pSet.get<int>("debug.trackingVerbosityLevel",0)),
    _smvlevel(pSet.get<int>("debug.steppingVerbosityLevel",0)),
    _checkFieldMap(pSet.get<int>("debug.checkFieldMap",0)),
    
    _generatorModuleLabel(pSet.get<string>("generatorModuleLabel", "")),
    _physVolHelper(),
<<<<<<< HEAD
    //_printPhysicsProcessSummary(pSet.get<bool>("debug.printPhysicsProcessSummary",false)),
=======
    _extMonFNALPixelSD(),
>>>>>>> First stage of changes to code to integrate G4MT into art3 - get rid of stashes
    _genEventBroker(_use_G4MT),

    _timer(make_unique<G4Timer>()),
    _realElapsed(0.),
    _systemElapsed(0.),
    _userElapsed(0.),
    _standardMu2eDetector((art::ServiceHandle<GeometryService>())->isStandardMu2eDetector()),
    _masterThreadIndex(_use_G4MT ? _nThreads : 0),
    _engine{art::ServiceHandle<SeedService>{}->getSeed()}
    {
        if((_generatorModuleLabel == art::InputTag()) && multiStagePars_.genInputHits().empty()) {
            throw cet::exception("CONFIG")
            << "Error: both generatorModuleLabel and genInputHits are empty - nothing to do!\n";
        }
            
        //get the right type of RunManager depending on whether we are in MT mode or not
        if (_use_G4MT) {
            _runManager.reset(new Mu2eG4MTRunManager());
        }
        else {
            _runManager.reset(new G4RunManager());
        }

    //In sequential mode, we need 1 SDHelper and it lives at index 0.
    //
    //In MT, mode we need one SDHelper for each Worker thread, plus one extra for the Master
    //in the ActionInitialization, each worker thread is given one of these SDHs to hold its SD-related data
    //the "0th" worker thread gets the "0th" element of the vector, etc
    //we give the "_nThreads" element to the Master thread through Mu2eG4World to setup the InstanceMap in the ctor of the SDH class
    //we need only one of these SDHs to declare to art the list of products that will be produced
    int nSDHelpersNeeded = ( _use_G4MT ) ? _nThreads+1 : 1;
    SensitiveDetectorHelpers.reserve(nSDHelpersNeeded);

    auto sd_pSet = pSet.get<fhicl::ParameterSet>("SDConfig", fhicl::ParameterSet());;
    for (int i = 0; i != nSDHelpersNeeded; ++i) {
      SensitiveDetectorHelpers.emplace_back(sd_pSet);
    }

    auto& collector = producesCollector();
    SensitiveDetectorHelpers.at(_masterThreadIndex).declareProducts(collector);

    produces<StatusG4>();
    produces<SimParticleCollection>();
     
    //******* This statement requires that the external libraries the module uses are thread-safe,
    //******* and that the data member members are used in a thread-safe manner
    async<art::InEvent>();

 
    //can we simplify this and directly declare the relevent products
    //rather than contructing these unneccesary object?
<<<<<<< HEAD
    stackingCuts_->declareProducts(collector);
    steppingCuts_->declareProducts(collector);
    commonCuts_->declareProducts(collector);
=======
//    stackingCuts_->declareProducts(this);
//    steppingCuts_->declareProducts(this);
//    commonCuts_->declareProducts(this);
>>>>>>> First stage of changes to code to integrate G4MT into art3 - get rid of stashes

    // Declare which products this module will read.
   
  
    
    if (_generatorModuleLabel != invalid_tag) {
      consumes<GenParticleCollection>(_generatorModuleLabel);
    }
    

    // The string "G4Engine" is magic; see the docs for RandomNumberGenerator.
    //createEngine( art::ServiceHandle<SeedService>()->getSeed(), "G4Engine");
        
        

} // end G4:G4(fhicl::ParameterSet const& pSet);


    void Mu2eG4::beginRun( art::Run &run, art::ProcessingFrame const& procFrame){

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
    if ( ncalls == 1 ) {
<<<<<<< HEAD
      stackingCuts_->finishConstruction(_originInWorld);
      steppingCuts_->finishConstruction(_originInWorld);
      commonCuts_->finishConstruction  (_originInWorld);

        //can only be run in single-threaded mode
      if( _checkFieldMap>0 && !(_use_G4MT)) generateFieldMap(_originInWorld,_checkFieldMap);
=======
      // Since the cuts that are put into the event are owned by the individual threads,
      // I do not think we need to do this.  Taking it out for now.
      //stackingCuts_->finishConstruction(originInWorld);
      //steppingCuts_->finishConstruction(originInWorld);
      //commonCuts_->finishConstruction(originInWorld);

      //can only be run in single-threaded mode
      if( _checkFieldMap>0 && !(_use_G4MT)) generateFieldMap(originInWorld,_checkFieldMap);
>>>>>>> First stage of changes to code to integrate G4MT into art3 - get rid of stashes

      if ( _exportPDTStart ) exportG4PDT( "Start:" );//once per job
    }//if ( ncalls == 1 )
    
}


void Mu2eG4::initializeG4( GeometryService& geom, art::Run const& run ){

    //if running in MT mode, set number of threads.
    //need to downcast the ptr to the RunManager, which was defined as a ptr to G4RunManager
    if (_use_G4MT) {
        dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->SetNumberOfThreads(_nThreads);
    }

    if ( _rmvlevel > 0 ) {
      mf::LogInfo logInfo("GEOM");
      logInfo << "Initializing Geant4 for " << run.id()
              << " with verbosity " << _rmvlevel << endl;
      logInfo << " Configured simParticleNumberOffset = "<< multiStagePars_.simParticleNumberOffset() << endl;
    }


    // Create user actions and register them with G4.
    G4VUserDetectorConstruction* allMu2e;

<<<<<<< HEAD
    // as mentioned above, we give the last element to the Master thread to setup the InstanceMap in the ctor of the SDH class
    
    if (_standardMu2eDetector) {
      geom.addWorldG4(*GeomHandle<Mu2eHall>());

      allMu2e = 
     	(new WorldMaker<Mu2eWorld>(std::make_unique<Mu2eWorld>(pset_, &(SensitiveDetectorHelpers.at(_masterThreadIndex))  ),
     				   std::make_unique<ConstructMaterials>(pset_)) );
    
      _originInWorld = (GeomHandle<WorldG4>())->mu2eOriginInWorld();
    }
    else {
      allMu2e =
	(new WorldMaker<Mu2eStudyWorld>(std::make_unique<Mu2eStudyWorld>(pset_, &(SensitiveDetectorHelpers.at(_masterThreadIndex)) ),
					std::make_unique<ConstructMaterials>(pset_)) );

      // non-Mu2e detector: the system origin os set to (0.,0.,0.); do not use geometry service for that
      _originInWorld = G4ThreeVector(0.0,0.0,0.0);
    }


    // in the non Mu2e detector we are working in the system with the
    // origin set to 0.,0.,0. and do not use geometry service for that
    //    originInWorld = (!_standardMu2eDetector) ? G4ThreeVector(0.0,0.0,0.0) : (GeomHandle<WorldG4>())->mu2eOriginInWorld();
=======
    //as mentioned above, we give the last element to the Master thread to setup the InstanceMap in the ctor of the SDH class
    if (standardMu2eDetector_) {
        allMu2e =
          (new WorldMaker<Mu2eWorld>(make_unique<Mu2eWorld>(pset_, &(SensitiveDetectorHelpers.at(_masterThreadIndex))  ),
                                       make_unique<ConstructMaterials>(pset_)) );
    }
    else {
        allMu2e =
        (new WorldMaker<Mu2eStudyWorld>(make_unique<Mu2eStudyWorld>(pset_, &(SensitiveDetectorHelpers.at(_masterThreadIndex)) ),
                                            make_unique<ConstructMaterials>(pset_)) );
    }

    // in the non-Mu2e detector we are working in the system with the
    // origin set to (0.,0.,0.) and do not use geometry service for that
    originInWorld = (!standardMu2eDetector_) ? G4ThreeVector(0.0,0.0,0.0) : (GeomHandle<WorldG4>())->mu2eOriginInWorld();
>>>>>>> First stage of changes to code to integrate G4MT into art3 - get rid of stashes

    preG4InitializeTasks(pset_);

    _runManager->SetVerboseLevel(_rmvlevel);
    _runManager->SetUserInitialization(allMu2e);

    physicsList_ = physicsListDecider(pset_);
    physicsList_->SetVerboseLevel(_rmvlevel);

    G4ParticleHPManager::GetInstance()->SetVerboseLevel(_rmvlevel);
    G4HadronicProcessStore::Instance()->SetVerbose(_rmvlevel);

    _runManager->SetUserInitialization(physicsList_);

    //this is where the UserActions are instantiated
    ActionInitialization* actioninit = new ActionInitialization(pset_,
                                                                SensitiveDetectorHelpers,
                                                                &_genEventBroker, &_physVolHelper,
                                                                _use_G4MT, _nThreads, _originInWorld,
                                                                mu2elimits_,
                                                                multiStagePars_.simParticleNumberOffset()
                                                                );

    //in MT mode, this is where BuildForMaster is called for master thread
    // in sequential mode, this is where Build() is called for main thread
    _runManager->SetUserInitialization(actioninit);

    if ( _rmvlevel > 0 ) {
      //  GetRunManagerType() returns an enum named RMType to indicate
      //  what kind of RunManager it is. RMType is defined as {
      //  sequentialRM, masterRM, workerRM }

       G4cout << __func__
              << " Before Initialize(), G4RunManagerType: "
              << _runManager->GetRunManagerType() << G4endl;
    }
    
    // Initialize G4 for this run.  Geometry and Physics are initialized.
    // IN MT mode, BeamOn(0) is called to do additional setup such as setting up the ScoringWorlds.
    _runManager->Initialize();

    if ( _rmvlevel > 0 ) {
       G4cout << __func__
              << " After Initialize(), G4RunManagerType:  "
              << _runManager->GetRunManagerType() << G4endl;
    }

} // end G4::initializeG4


    void Mu2eG4::beginSubRun(art::SubRun& sr, art::ProcessingFrame const& procFrame)
{

    
}


// Create one G4 event and copy its output to the art::event.
    void Mu2eG4::produce(art::Event& event, art::ProcessingFrame const& procFrame) {

    if (_use_G4MT && event.id().event() == 1) cout << "\n*-*-*-*-*-*- You are running "
       << dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->GetNumberOfThreads() << " threads. -*-*-*-*-*-*\n" << endl;


    
    //******** these are per art::event quantities ********

    // Handle to the generated particles; need when building art::Ptr to a GenParticle.
    art::Handle<GenParticleCollection> gensHandle;
    if(!(_generatorModuleLabel == art::InputTag())) {
        event.getByLabel(_generatorModuleLabel, gensHandle);
    }
    
    // StepPointMCCollection of input hits from the previous simulation stage
    HitHandles genInputHits;
    for(const auto& i : multiStagePars_.genInputHits()) {
        genInputHits.emplace_back(event.getValidHandle<StepPointMCCollection>(i));
    }

    // WILL LIKLEY NEED THIS TO GIVE TO THE PEREVENTOBJECTSMANAGER VIA the GENEVENTBROKER?
    // ProductID for the SimParticleCollection.
    art::ProductID simPartId(event.getProductID<SimParticleCollection>());
    art::EDProductGetter const* simProductGetter = event.productGetter(simPartId);
    
    // Create empty data products.
    unique_ptr<SimParticleCollection>      simParticles(      new SimParticleCollection);
    unique_ptr<ExtMonFNALSimHitCollection> extMonFNALHits(    new ExtMonFNALSimHitCollection);
    
    
    _genEventBroker.loadEvent(genInputHits, simPartId, &event, _generatorModuleLabel, simProductGetter);


    //TEMP FIX, we need to get the number of events to process from fcl file
    int number_of_events_to_process = 100;

    // Run G4 for this event and access the completed event.
    BeamOnDoOneArtEvent( event.id().event(), number_of_events_to_process );

    _genEventBroker.setEventPtrToZero();

/////////////////////////////////////////////////////////////////////////////////////
    
    if (simParticles == nullptr) {
        numExcludedEvents++;
    }
    else {
        
        //THESE ARE CURRENTLY IN THE EVENTACTION
        //event.put(move(g4stat)); THIS data product doesn't exist in this context.
        //event.put(move(simParticles));
        
        //need something similar for the SensitiveDetector hits
        
        
//TAKING THESE OUT TEMPORARILY
//        stackingCuts_->put(event);
//        steppingCuts_->put(event);
//        commonCuts_->put(event);
    
    }
    
}//end Mu2eG4::produce


// Tell G4 that this run is over.
    void Mu2eG4::endRun(art::Run & run, art::ProcessingFrame const& procFrame){

        BeamOnEndRun();

        G4cout << "at endRun: numExcludedEvents = " << numExcludedEvents << G4endl;

}


    void Mu2eG4::endJob(art::ProcessingFrame const& procFrame){

    if ( _exportPDTEnd ) exportG4PDT( "End:" );

    // Yes, these are named endRun, but they are really endJob actions.
    _physVolHelper.endRun();
}





/**************************************************************
                    FUNCTION DEFINITIONS
 **************************************************************/

// Do the "begin run" parts of BeamOn.
void Mu2eG4::BeamOnBeginRun( unsigned int runNumber){

  if ( _rmvlevel > 0 ) {
    G4cout << __func__ << ": runNumber: " << runNumber << G4endl;
  }


        _runManager->SetRunIDCounter(runNumber);

        bool cond = _runManager->ConfirmBeamOnCondition();
        if(!cond){
            // throw here
            return;
        }

        _realElapsed   = 0.;
        _systemElapsed = 0.;
        _userElapsed   = 0.;

        //NOTE: the data member G4RunManager::numberOfEventToBeProcessed would have been set from within the G4 call to BeamOn()
        //BEFORE the call the ConstructScoringWorlds() with this call:
        //_runManager->SetNumberOfEventsToBeProcessed(numberOfEventsToBeProcessed);
        //Since we have deconstructed the BeamOn() function, we must set this data member ourselves.
        //In the pre-MT version of this code, this call was made here, immediately before the call _runManager->ConstructScoringWorlds();
        //It has been moved, because at 'BeginRun' time, the value of this variable cannot be known without using a FHiCL parameter.
        //It was decided not to define it though a FHiCL parameter, but rather within the code.
        //
        //In SEQUENTIAL mode, this variable doesn't seem to serve any purpose other than to be used in the call
        //currentRun->SetNumberOfEventToBeProcessed(numberOfEventToBeProcessed); in G4RunManager::RunInitialization().
        //Furthermore, at least up through G4 v10.4.02 the data member G4Run::numberOfEventToBeProcessed is only informative.
        //It is not needed anywhere in the G4 code.  So, it is OK to not set this variable at this point.
        //
        //For MT mode, G4MTRunManager::numberOfEventToBeProcessed is used to distribute events to the threads, among other things,
        //and must have the value of EITHER the size of the GenParticle Stash when using single-threaded art,
        //OR the total number of events to be simulated when using multithreaded art.
        //It must be set before the call to G4MTRunManager::InitializeEventLoop.
        //We can get the size of the GenParticle Stash from WITHIN the code with the call to _genEventBroker.getStashSize(),
        //but this cannot be called until after _genEventBroker.loadEvent() is called.
        //So, we get this variable from the GenParticle Stash after the call to _genEventBroker.loadEvent().  We initialize the
        //EventStash to the correct size using it, and we pass it into G4 in the call BeamOnDoOneArtEvent().

        _runManager->ConstructScoringWorlds();
        _runManager->RunInitialization();

}//BeamOnBeginRun


// Do the "per event" part of DoEventLoop.
void Mu2eG4::BeamOnDoOneArtEvent( int eventNumber, G4int num_events, const char* macroFile, G4int n_select){

  if ( _rmvlevel > 1 ) {
    G4cout << __func__ << ": eventNumber: " << eventNumber << ", num_events: " << num_events << G4endl;
  }



        if (_use_G4MT)//MT mode
        {
            //this is where the events are actually processed
            //num_events is # of G4 events processed per art event

            //NOTE: G4MTRunManager::WaitForEndEventLoopWorkers() is a protected function in G4 code, so we MUST have our own MTRunManager
            //in order to access this function

            dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->SetNumberOfEventsToBeProcessed(num_events);
            dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->InitializeEventLoop(num_events,macroFile,n_select);
            dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->Mu2eG4WaitForEndEventLoopWorkers();
            //dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->Mu2eG4TerminateWorkers();//DID NOT WORK BEFORE

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

}//BeamOnDoOneArtEvent


// Do the "end of run" parts of DoEventLoop and BeamOn.
void Mu2eG4::BeamOnEndRun(){


    if (_use_G4MT)//MT mode
    {
        dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->Mu2eG4RunTermination();
    }
    else//sequential mode
    {
        _runManager->TerminateEventLoop();
        _runManager->RunTermination();

 
    }
}//BeamOnEndRun


} // End of namespace mu2e

DEFINE_ART_MODULE(mu2e::Mu2eG4);
