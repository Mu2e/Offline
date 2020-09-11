//
// Check the output of FilterG4Out
//
// $Id: CompressionCheck0_module.cc,v 1.2 2013/10/21 20:44:04 genser Exp $
// $Author: genser $
// $Date: 2013/10/21 20:44:04 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <string>
#include <set>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

// Mu2e includes.
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"
#include "GeneralUtilities/inc/Binning.hh"

// Root includes.
#include "TH1F.h"

using namespace std;

namespace{

  TH1F* bookTH1F( art::TFileService& tfs,
                    std::string const& name,
                    std::string const& title,
                    ::Binning const& bins){
    TH1F* val = tfs.make<TH1F>( name.c_str(), title.c_str(), bins.nbins(), bins.low(), bins.high() );
    return val;
  }
}

namespace mu2e {

  class CompressionCheck0 : public art::EDAnalyzer {
  public:

    explicit CompressionCheck0(fhicl::ParameterSet const& pset);
    virtual ~CompressionCheck0() { }

    // The framework calls this at the start of the job.
    virtual void beginJob();

    // The framework calls this for each event.
    void analyze(const art::Event& e);

  private:

    // Module label of the g4 module that made the hits.
    std::string g4ModuleLabel_;
    std::string filterModuleLabel_;

    std::string trackerStepsName_;
    std::string caloStepsName_;
    std::string croStepsName_;
    std::string vdStepsName_;

    art::InputTag trackerTag1_;
    art::InputTag trackerTag2_;

    art::InputTag caloTag1_;
    art::InputTag caloTag2_;

    art::InputTag croTag1_;
    art::InputTag croTag2_;

    art::InputTag vdTag1_;
    art::InputTag vdTag2_;

    art::InputTag simsTag1_;
    art::InputTag simsTag2_;

    // Bin definitions: must be specified in the 
    Binning trackerStepBins_;
    Binning caloStepBins_;
    Binning croStepBins_;
    Binning vdStepBins_;
    Binning simsBeforeBins_;
    Binning simsAfterBins_;

    // Pointers to histograms, ntuples.
    TH1F* hTrackerSteps1_;
    TH1F* hTrackerSteps2_;
    TH1F* hTrackerStepsD_;
    TH1F* hCaloSteps1_;
    TH1F* hCaloSteps2_;
    TH1F* hCaloStepsD_;
    TH1F* hcroSteps1_;
    TH1F* hcroSteps2_;
    TH1F* hcroStepsD_;
    TH1F* hvdSteps1_;
    TH1F* hvdSteps2_;
    TH1F* hnSims1_;
    TH1F* hnSims2_;

  };

  CompressionCheck0::CompressionCheck0(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    g4ModuleLabel_    (pset.get<string>("g4ModuleLabel")),
    filterModuleLabel_(pset.get<string>("filterModuleLabel")),

    trackerStepsName_( StepInstanceName(StepInstanceName::tracker).name() ),
    caloStepsName_   ( StepInstanceName(StepInstanceName::calorimeter).name() ),
    croStepsName_    ( StepInstanceName(StepInstanceName::calorimeterRO).name() ),
    vdStepsName_     ( StepInstanceName(StepInstanceName::virtualdetector).name() ),

    trackerTag1_( g4ModuleLabel_     + ":" + trackerStepsName_ ),
    trackerTag2_( filterModuleLabel_ + ":" + trackerStepsName_ ),

    caloTag1_( g4ModuleLabel_     + ":" + caloStepsName_ ),
    caloTag2_( filterModuleLabel_ + ":" + caloStepsName_ ),

    croTag1_( g4ModuleLabel_     + ":" + croStepsName_ ),
    croTag2_( filterModuleLabel_ + ":" + croStepsName_ ),

    vdTag1_( g4ModuleLabel_     + ":" + vdStepsName_ ),
    vdTag2_( filterModuleLabel_ + ":" + vdStepsName_ ),

    simsTag1_( g4ModuleLabel_              ),
    simsTag2_( filterModuleLabel_  + ":s0" ),

    trackerStepBins_(pset.get<Binning>("trackerStepBins")),
    caloStepBins_   (pset.get<Binning>("caloStepBins")),
    croStepBins_    (pset.get<Binning>("croStepBins")),
    vdStepBins_     (pset.get<Binning>("vdStepBins")),
    simsBeforeBins_  (pset.get<Binning>("simsBeforeBins")),
    simsAfterBins_   (pset.get<Binning>("simsAfterBins")),

    // Histograms
    hTrackerSteps1_(0),
    hTrackerSteps2_(0),
    hTrackerStepsD_(0),
    hCaloSteps1_(0),
    hCaloSteps2_(0),
    hCaloStepsD_(0),
    hcroSteps1_(0),
    hcroSteps2_(0),
    hcroStepsD_(0),
    hvdSteps1_(0),
    hvdSteps2_(0),

    hnSims1_(0),
    hnSims2_(0){

  }

  // At the start of the job, book histograms.
  void CompressionCheck0::beginJob(){

    // Get a handle to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    hTrackerSteps1_ = bookTH1F( *tfs, "trackerSteps1", "Tracker Steps Before", trackerStepBins_ );
    hTrackerSteps2_ = bookTH1F( *tfs, "trackerSteps2", "Tracker Steps After", trackerStepBins_ );
    hTrackerStepsD_ = tfs->make<TH1F>( "trackerStepsD", "Tracker Steps Delta",  3, -1.,   1. );

    hCaloSteps1_ = bookTH1F( *tfs, "caloSteps1", "Calo Steps Before",  caloStepBins_ );
    hCaloSteps2_ = bookTH1F( *tfs, "caloSteps2", "Calo Steps After",  caloStepBins_ );
    hCaloStepsD_ = tfs->make<TH1F>( "caloStepsD", "Calo Steps Delta",    3, -1.,   1. );

    hcroSteps1_ = bookTH1F( *tfs, "croSteps1", "Calo Readout Steps Before", croStepBins_ );
    hcroSteps2_ = bookTH1F( *tfs, "croSteps2", "Calo Readout Steps After", croStepBins_ );
    hcroStepsD_ = tfs->make<TH1F>( "croStepsD", "Calo Readout Steps Delta",    3, -1.,   1. );

    hvdSteps1_ = bookTH1F( *tfs, "vdSteps1", "vd Steps Before",  vdStepBins_ );
    hvdSteps2_ = bookTH1F( *tfs, "vdSteps2", "vd Steps After",   vdStepBins_ );

    hnSims1_ = bookTH1F( *tfs, "nSims1", "Number of SimParticles Before", simsBeforeBins_ );
    hnSims2_ = bookTH1F( *tfs, "nSims2", "Number of SimParticles After",  simsAfterBins_  );

  } // end beginJob

  // For each event, look at tracker hits and calorimeter hits.
  void CompressionCheck0::analyze(const art::Event& event) {

    // Tracker
    art::Handle<StepPointMCCollection> tsteps1;
    event.getByLabel(trackerTag1_, tsteps1);
    auto tsteps2 = event.getValidHandle<StepPointMCCollection>(trackerTag2_);

    if ( tsteps1.isValid() ) {
      hTrackerSteps1_->Fill( tsteps1->size() );
      hTrackerStepsD_->Fill( tsteps2->size()-tsteps1->size() );
    }
    hTrackerSteps2_->Fill( tsteps2->size() );

    // Calorimeter
    art::Handle<StepPointMCCollection> csteps1;
    event.getByLabel(caloTag1_, csteps1);
    auto csteps2 = event.getValidHandle<StepPointMCCollection>(caloTag2_);

    if ( csteps1.isValid() ) {
      hCaloSteps1_->Fill( csteps1->size() );
      hCaloStepsD_->Fill( csteps2->size()-csteps1->size() );
    }
    hCaloSteps2_->Fill( csteps2->size() );

    // Calorimeter readout
    art::Handle<StepPointMCCollection> rsteps1;
    event.getByLabel(croTag1_, rsteps1);
    auto rsteps2 = event.getValidHandle<StepPointMCCollection>(croTag2_);

    if ( rsteps1.isValid() ) {
      hcroSteps1_->Fill( rsteps1->size() );
      hcroStepsD_->Fill( rsteps2->size()-rsteps1->size() );
    }
    hcroSteps2_->Fill( rsteps2->size() );

    // Virtual detector
    art::Handle<StepPointMCCollection> vsteps1;
    event.getByLabel(vdTag1_, vsteps1);
    auto vsteps2 = event.getValidHandle<StepPointMCCollection>(vdTag2_);

    if ( vsteps1.isValid() ) {
      hvdSteps1_->Fill( vsteps1->size() );
    }
    hvdSteps2_->Fill( vsteps2->size() );

    // SimParticles
    art::Handle<SimParticleCollection> sims1;
    event.getByLabel( simsTag1_, sims1);
    auto sims2 = event.getValidHandle<SimParticleCollection>( simsTag2_);

    if ( sims1.isValid() ){
      hnSims1_->Fill(sims1->size());
    }
    hnSims2_->Fill(sims2->size());

  }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CompressionCheck0);
