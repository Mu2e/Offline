//
// Plugin to show how to use the SimParticlesWithHits class.
//
// $Id: Summary01_module.cc,v 1.10 2013/10/21 20:44:04 genser Exp $
// $Author: genser $
// $Date: 2013/10/21 20:44:04 $
//
// Original author Rob Kutschke.
//

// C++ includes.
#include <iostream>
#include <string>
#include <map>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

// Mu2e includes.
#include "TrackerGeom/inc/Tracker.hh"
#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"

// ROOT includes
#include "TH1F.h"

using namespace std;

namespace mu2e {

  class Summary01 : public art::EDAnalyzer {
  public:
    explicit Summary01(fhicl::ParameterSet const& pset):
      art::EDAnalyzer(pset),
      g4ModuleLabel_(pset.get<std::string>("g4ModuleLabel")),
      strawHitMakerModuleLabel_(pset.get<std::string>("strawHitMakerModuleLabel")),
      trackerStepPoints_(pset.get<std::string>("trackerStepPoints")),
      calorimeterStepPoints_(pset.get<std::string>("calorimeterStepPoints")),
      calorimeterROStepPoints_(pset.get<std::string>("calorimeterROStepPoints")),
      crvStepPoints_(pset.get<std::string>("crvStepPoints")),
      minEnergyDep_(pset.get<double>("minEnergyDep")),
      minHits_(pset.get<unsigned>("minHits")),
      simsPlotMax_(pset.get<double>("simsPlotMax",400.)),
      productionMode_(pset.get<bool>("productionMode",true)),
      deltaRayParentId_(),
      deltaRayVolumeId_(),
      nBadG4_(0),

      // Histograms
      hG4Status_(0),
      hDeltaRayEKine0_(0),
      hDeltaRayEKine1_(0),
      hDeltaRayEKine2_(0),
      hSimsSize_(0),
      hTrkSimsSize_(0),
      hCalSimsSize_(0),
      hTrkSimsRatio_(0),
      hCalSimsRatio_(0){
    }
    virtual ~Summary01() { }

    void beginJob();
    void endJob();

    void analyze( art::Event const& e );

    void endRun( art::Run const& );


  private:

    // Label of the modules that created the data products.
    std::string g4ModuleLabel_;
    std::string strawHitMakerModuleLabel_;

    // Names of the interesting StepPoint collections.
    std::string trackerStepPoints_;
    std::string calorimeterStepPoints_;
    std::string calorimeterROStepPoints_;
    std::string crvStepPoints_;

    // Cuts used inside SimParticleWithHits:
    //  - drop hits with too little energy deposited.
    //  - drop SimParticles with too few hits.
    double minEnergyDep_;
    size_t minHits_;

    // Maximum size of the plots to study compression of the SimParticleCollection.
    double simsPlotMax_;

    // If true then development code is turned off and verbosity is dropped.
    bool productionMode_;

    std::map<int,int> deltaRayParentId_;
    std::map<int,int> deltaRayVolumeId_;

    int nBadG4_;

    TH1F* hG4Status_;

    TH1F* hDeltaRayEKine0_;
    TH1F* hDeltaRayEKine1_;
    TH1F* hDeltaRayEKine2_;

    TH1F* hSimsSize_;
    TH1F* hTrkSimsSize_;
    TH1F* hCalSimsSize_;
    TH1F* hTrkSimsRatio_;
    TH1F* hCalSimsRatio_;

    void deltaRaySpectra( art::Event const& );
    void compressSims( art::Event const& );

  };

  void
  Summary01::beginJob( ){

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    hG4Status_ = tfs->make<TH1F>( "hG4Status", "Non-zero Values of G4 Status word", 40, 0.,  40.   );

    hDeltaRayEKine0_ = tfs->make<TH1F>( "hDeltaRayEKine0", "Kinetic Energy of delta-ray;(MeV)", 100, 0.,  100.   );
    hDeltaRayEKine1_ = tfs->make<TH1F>( "hDeltaRayEKine1", "Kinetic Energy of delta-ray;(MeV)", 100, 0.,    5.   );
    hDeltaRayEKine2_ = tfs->make<TH1F>( "hDeltaRayEKine2", "Kinetic Energy of delta-ray;(MeV)", 100, 0.,    0.05 );


    hSimsSize_     = tfs->make<TH1F>( "hSimsSize",     "Size of SimParticleColleciton",            100, 0.,  simsPlotMax_ );
    hTrkSimsSize_  = tfs->make<TH1F>( "hTrkSimsSize",  "SimParticleColleciton: size for non-Cal",  100, 0.,  simsPlotMax_ );
    hCalSimsSize_  = tfs->make<TH1F>( "hCalSimsSize",  "SimParticleColleciton: size for all Hits", 100, 0.,  simsPlotMax_ );

    // Make sure that 100% is in the plot.
    hTrkSimsRatio_ = tfs->make<TH1F>( "hTrkSimsRatio", "SimParticleColleciton: Ratio for non-Cal",  101, 0.,  1.01 );
    hCalSimsRatio_ = tfs->make<TH1F>( "hCalSimsRatio", "SimParticleColleciton: Ratio for all Hits", 101, 0.,  1.01 );

  }

  void
  Summary01::analyze(art::Event const& event ) {

    // Skip event if G4 did not complete OK.
    art::Handle<StatusG4> g4StatusHandle;
    event.getByLabel( g4ModuleLabel_, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;
    if ( g4Status.status() > 1 ){
      ++nBadG4_;
      hG4Status_->Fill( g4Status.status() );
      mf::LogError("G4")
        << "Summary01::analyze skipping event because of bad G4 status.\n"
        << g4Status;
      return;
    }

    // Ask the event to give us a "handle" to the requested hits.
    //art::Handle<StepPointMCCollection> stepsHandle;
    //event.getByLabel( g4ModuleLabel_, trackerStepPoints_,stepsHandle);
    //StepPointMCCollection const& steps = *stepsHandle;

    //art::Handle<SimParticleCollection> simsHandle;
    //event.getByLabel( g4ModuleLabel_, simsHandle);
    //SimParticleCollection const& sims = *simsHandle;

    art::Handle<PhysicalVolumeInfoCollection> volsHandle;
    event.getRun().getByLabel( g4ModuleLabel_, volsHandle);
    PhysicalVolumeInfoCollection const& vols = *volsHandle;

    // Fill plots of the spectra of delta rays.
    deltaRaySpectra( event );
    compressSims( event );

    // Stuff below here is under development.
    if ( productionMode_ ) return;

    // Construct an object that ties together all of the simulated particle and hit info.
    SimParticlesWithHits simsWithHits( event,
                                       g4ModuleLabel_,
                                       strawHitMakerModuleLabel_,
                                       trackerStepPoints_,
                                       minEnergyDep_,
                                       minHits_ );

    typedef SimParticlesWithHits::map_type map_type;

    for ( map_type::const_iterator i=simsWithHits.begin();
          i != simsWithHits.end(); ++i ){

      // All information about this SimParticle
      SimParticleInfo const& simInfo = i->second;

      // Information about StrawHits that belong on this SimParticle.
      vector<StrawHitMCInfo> const& infos = simInfo.strawHitInfos();

      if ( infos.size() > 5 && i->first != SimParticleCollection::key_type(1) ){

        SimParticle const& sim = simInfo.simParticle();

        cout << "SimParticle: "
             << " Event: " << event.id().event()
             << " Track: " << i->first
             << " PdgId: " << sim.pdgId()
             << " |p|: "   << sim.startMomentum().vect().mag()
             << " Hits: "  << infos.size()
             << " | "
             << vols[sim.startVolumeIndex()] << " "
             << vols[sim.endVolumeIndex()]
             << endl;
      }

    }

  } // end of ::analyze.

  void
  Summary01::deltaRaySpectra( art::Event const& event ){

    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel( g4ModuleLabel_, simsHandle);
    SimParticleCollection const& sims = *simsHandle;

    // Plot the energy spectrum of delta rays.
    for ( SimParticleCollection::const_iterator i=sims.begin();
          i!=sims.end(); ++i ){
      SimParticle const& sim = i->second;

      if ( sim.pdgId()        != PDGCode::e_minus   ) continue;
      if ( !sim.madeInG4()                          ) continue;
      if ( sim.creationCode() != ProcessCode::eIoni ) continue;

      SimParticle const& parent = sims[sim.parentId()];
      if ( parent.pdgId() == PDGCode::gamma ) continue;

      // Kinetic energy
      double ek = sim.startMomentum().e()-sim.startMomentum().restMass();

      // Fill histograms
      hDeltaRayEKine0_->Fill(ek);
      hDeltaRayEKine1_->Fill(ek);
      hDeltaRayEKine2_->Fill(ek);

      // Count how many delta rays came from each parent.
      ++deltaRayParentId_[parent.pdgId()];

      ++deltaRayVolumeId_[sim.startVolumeIndex()];
    }

  } // end Summary01::deltaRaySpectra


  // Out of class utility function used in compressSims.
  // Navigate from each StepPointMC back to the PrimaryGeneratorAction particle
  // and mark each track along the way as being on the path to a hit.
  void markChainToHit( map<SimParticleCollection::key_type,int>& used,
                       StepPointMCCollection const& steps,
                       SimParticleCollection const& sims ){

    typedef SimParticleCollection::key_type key_type;

    for ( size_t i=0; i<steps.size(); ++i ){
      StepPointMC const& step = steps.at(i);
      key_type id = step.trackId();

      // Mark this track has being in the chain to a hit.
      ++used[id];

      // Navigate back to the generated particle
      SimParticle const* sim = &sims[id];
      while ( sim->hasParent() ){
        key_type parentKey = sim->parentId();
        ++used[parentKey];
        sim = &sims[parentKey];
      }
    }
  }


  // Investigate how much space we can save by compressing the SimParticleCollection.
  void
  Summary01::compressSims( art::Event const& event ){

    art::Handle<StepPointMCCollection> trkStepsHandle;
    event.getByLabel( g4ModuleLabel_, trackerStepPoints_,trkStepsHandle);
    StepPointMCCollection const& trkSteps = *trkStepsHandle;

    art::Handle<StepPointMCCollection> calStepsHandle;
    event.getByLabel( g4ModuleLabel_, calorimeterStepPoints_,calStepsHandle);
    StepPointMCCollection const& calSteps = *calStepsHandle;

    art::Handle<StepPointMCCollection> cROStepsHandle;
    event.getByLabel( g4ModuleLabel_, calorimeterROStepPoints_,cROStepsHandle);
    StepPointMCCollection const& cROSteps = *cROStepsHandle;

    art::Handle<StepPointMCCollection> crvStepsHandle;
    event.getByLabel( g4ModuleLabel_, calorimeterROStepPoints_,crvStepsHandle);
    StepPointMCCollection const& crvSteps = *crvStepsHandle;

    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel( g4ModuleLabel_, simsHandle);
    SimParticleCollection const& sims = *simsHandle;

    // Count how often each track is on the path from
    // a StepPointMC to a generated particle.
    typedef SimParticleCollection::key_type key_type;
    map<key_type,int> used;

    // Contribution from tracker StepPointMCs.
    markChainToHit( used, trkSteps, sims);
    //size_t size1 = used.size();

    // Contribution from CRV StepPointMCs.
    markChainToHit( used, crvSteps, sims);
    //size_t size2 = used.size();

    // Contribution from calorimeter RO StepPointMCs.
    markChainToHit( used, cROSteps, sims);
    size_t size3 = used.size();

    // Contribution from calorimeter crystal StepPointMCs.
    // This is the big one.
    markChainToHit( used, calSteps, sims);
    size_t size4 = used.size();

    // Fill histograms to document this.
    size_t simsSize=sims.size();
    hSimsSize_->Fill(simsSize);
    hTrkSimsSize_->Fill(size3);
    hCalSimsSize_->Fill(size4);
    hTrkSimsRatio_->Fill( double(size3)/double(simsSize));
    hCalSimsRatio_->Fill( double(size4)/double(simsSize));

  }

  void Summary01::endRun( art::Run const& run ){

    // Skip immature verbose printout.
    if ( productionMode_ ) return;

    art::Handle<PhysicalVolumeInfoCollection> volsHandle;
    run.getByLabel( g4ModuleLabel_, volsHandle);
    PhysicalVolumeInfoCollection const& vols = *volsHandle;

    cout << "\n Creation volumes of Delta Rays: " << endl;
    for ( std::map<int,int>::const_iterator i=deltaRayVolumeId_.begin();
          i !=deltaRayVolumeId_.end(); ++i ){
      cout << " Volume: "
           << i->first << " "
           << vols[i->first] << " "
           << i->second
           << endl;
    }

    cout << "\n Parents of Delta Rays: " << endl;
    for ( std::map<int,int>::const_iterator i=deltaRayParentId_.begin();
          i !=deltaRayParentId_.end(); ++i ){
      cout << " "
           << i->first << " "
           << i->second
           << endl;
    }

  } // Summary01::endRun


  void Summary01::endJob(){
    if ( nBadG4_ > 0 ){
      cout << "\nNumber of events with failed G4: " << nBadG4_ << endl;
    } else{
      cout << "\nAll events completed G4 correctly" << endl;
    }
  } // end Summary01::endJob

}

using mu2e::Summary01;
DEFINE_ART_MODULE(Summary01);
