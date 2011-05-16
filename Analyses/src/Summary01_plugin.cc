//
// Plugin to show how to use the SimParticlesWithHits class.
//
// $Id: Summary01_plugin.cc,v 1.1 2011/05/16 00:23:27 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/16 00:23:27 $
//
// Original author Rob Kutschke.
//

// C++ includes.
#include <iostream>
#include <string>
#include <map>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Services/interface/TFileService.h"

// Mu2e includes.
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "ToyDP/inc/StatusG4.hh"

// ROOT includes
#include "TH1F.h"

using namespace std;

namespace mu2e {

  class Summary01 : public edm::EDAnalyzer {
  public:
    explicit Summary01(edm::ParameterSet const& pset):
      g4ModuleLabel_(pset.getParameter<std::string>("g4ModuleLabel")),
      strawHitMakerModuleLabel_(pset.getParameter<std::string>("strawHitMakerModuleLabel")),
      trackerStepPoints_(pset.getParameter<std::string>("trackerStepPoints")),
      calorimeterStepPoints_(pset.getParameter<std::string>("calorimeterStepPoints")),
      calorimeterROStepPoints_(pset.getParameter<std::string>("calorimeterROStepPoints")),
      crvStepPoints_(pset.getParameter<std::string>("crvStepPoints")),
      minEnergyDep_(pset.getParameter<double>("minEnergyDep")),
      minHits_(pset.getParameter<uint32_t>("minHits")),
      simsPlotMax_(pset.getUntrackedParameter<double>("simsPlotMax",400.)),
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

    void beginJob(edm::EventSetup const& );
    void endJob();

    void analyze( edm::Event const& e, edm::EventSetup const&);

    void endRun( edm::Run const&, edm::EventSetup const&);


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

    void deltaRaySpectra( edm::Event const& );
    void compressSims( edm::Event const& );

  };

  void 
  Summary01::beginJob(edm::EventSetup const& ){

    // Get access to the TFile service.
    edm::Service<edm::TFileService> tfs;

    hG4Status_ = tfs->make<TH1F>( "hG4Status", "Non-zero Values of G4 Status word", 40, 0.,  40.   );

    hDeltaRayEKine0_ = tfs->make<TH1F>( "hDeltaRayEKine0", "Kinetic Energy of delta-ray;(MeV)", 100, 0.,  100.   );
    hDeltaRayEKine1_ = tfs->make<TH1F>( "hDeltaRayEKine1", "Kinetic Energy of delta-ray;(MeV)", 100, 0.,    5.   );
    hDeltaRayEKine2_ = tfs->make<TH1F>( "hDeltaRayEKine2", "Kinetic Energy of delta-ray;(MeV)", 100, 0.,    0.05 );


    hSimsSize_     = tfs->make<TH1F>( "hSimsSize",     "Size of SimParticleColleciton",            100, 0.,  simsPlotMax_ );
    hTrkSimsSize_  = tfs->make<TH1F>( "hTrkSimsSize",  "SimParticleColleciton: size for non-Cal",  100, 0.,  simsPlotMax_ );
    hCalSimsSize_  = tfs->make<TH1F>( "hCalSimsSize",  "SimParticleColleciton: size for all Hits", 100, 0.,  simsPlotMax_ );
    hTrkSimsRatio_ = tfs->make<TH1F>( "hTrkSimsRatio", "SimParticleColleciton: Ratio for non-Cal",  100, 0.,  1.  );
    hCalSimsRatio_ = tfs->make<TH1F>( "hCalSimsRatio", "SimParticleColleciton: Ratio for all Hits", 100, 0.,  1.  );

  }

  void
  Summary01::analyze(edm::Event const& event, edm::EventSetup const&) {

    // Skip event if G4 did not complete OK.
    edm::Handle<StatusG4> g4StatusHandle;
    event.getByLabel( g4ModuleLabel_, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;
    if ( g4Status.status() != 0 ){
      ++nBadG4_;
      hG4Status_->Fill( g4Status.status() );
      edm::LogError("G4")
        << "Summary01::analyze skipping event because of bad G4 status\n"
        << g4Status;
      return;
    }

    // Ask the event to give us a "handle" to the requested hits.
    edm::Handle<StepPointMCCollection> stepsHandle;
    event.getByLabel( g4ModuleLabel_, trackerStepPoints_,stepsHandle);
    StepPointMCCollection const& steps = *stepsHandle;

    edm::Handle<SimParticleCollection> simsHandle;
    event.getByLabel( g4ModuleLabel_, simsHandle);
    SimParticleCollection const& sims = *simsHandle;

    edm::Handle<PhysicalVolumeInfoCollection> volsHandle;
    event.getRun().getByLabel( g4ModuleLabel_, volsHandle);
    PhysicalVolumeInfoCollection const& vols = *volsHandle;

    // Fill plots of the spectra of delta rays.
    deltaRaySpectra( event );
    compressSims( event );

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
  Summary01::deltaRaySpectra( edm::Event const& event ){

    edm::Handle<SimParticleCollection> simsHandle;
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
  Summary01::compressSims( edm::Event const& event ){

    edm::Handle<StepPointMCCollection> trkStepsHandle;
    event.getByLabel( g4ModuleLabel_, trackerStepPoints_,trkStepsHandle);
    StepPointMCCollection const& trkSteps = *trkStepsHandle;

    edm::Handle<StepPointMCCollection> calStepsHandle;
    event.getByLabel( g4ModuleLabel_, calorimeterStepPoints_,calStepsHandle);
    StepPointMCCollection const& calSteps = *calStepsHandle;

    edm::Handle<StepPointMCCollection> cROStepsHandle;
    event.getByLabel( g4ModuleLabel_, calorimeterROStepPoints_,cROStepsHandle);
    StepPointMCCollection const& cROSteps = *cROStepsHandle;

    edm::Handle<StepPointMCCollection> crvStepsHandle;
    event.getByLabel( g4ModuleLabel_, calorimeterROStepPoints_,crvStepsHandle);
    StepPointMCCollection const& crvSteps = *crvStepsHandle;

    edm::Handle<SimParticleCollection> simsHandle;
    event.getByLabel( g4ModuleLabel_, simsHandle);
    SimParticleCollection const& sims = *simsHandle;

    // Count how often each track is on the path from 
    // a StepPointMC to a generated particle.
    typedef SimParticleCollection::key_type key_type;
    map<key_type,int> used;

    // Contribution from tracker StepPointMCs.
    markChainToHit( used, trkSteps, sims);
    size_t size1 = used.size();

    // Contribution from CRV StepPointMCs.
    markChainToHit( used, crvSteps, sims);
    size_t size2 = used.size();

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

    cout << "Compare: " 
         << sims.size() << " " 
         << size1 << " "
         << size2 << " "
         << size3 << " "
         << size4 << " "
         << endl;

  }

  void Summary01::endRun( edm::Run const& run, edm::EventSetup const&){

    edm::Handle<PhysicalVolumeInfoCollection> volsHandle;
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
DEFINE_FWK_MODULE(Summary01);
