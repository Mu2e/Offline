//
// Select events with a minimum number of StepPointMC's in various detectors.
// $Id: MinimumHits_module.cc,v 1.2 2011/07/17 02:28:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/07/17 02:28:23 $
//
// Contact person Rob Kutschke.
//

// Mu2e includes.
#include "Analyses/inc/DiagnosticsG4.hh"
#include "Analyses/inc/GeneratorSummaryHistograms.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"

// Framework includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

// Root includes
#include "TH1F.h"

// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <iostream>

using namespace std;

namespace mu2e {

  class MinimumHits : public art::EDFilter {
  public:
    explicit MinimumHits(fhicl::ParameterSet const& pset);
    virtual ~MinimumHits() { }

    bool filter( art::Event& event);

    virtual bool beginRun(art::Run &run);

  private:

    // Module label of the g4 module that made the generated particles
    std::string generatorModuleLabel_;

    // Module label of the module that made the StepPointMCCollections.
    std::string g4ModuleLabel_;

    // Module labels of the modules that made the collection of reco hits.
    std::string strawHitMakerLabel_;
    std::string crystalHitMakerLabel_;

    // Instance names of the StepPointMCCollections.
    std::string trackerStepPoints_;
    std::string caloStepPoints_;
    std::string caloROStepPoints_;
    std::string foilStepPoints_;
    std::string crvStepPoints_;
    std::string vDetStepPoints_;

    // Histogram pointers.
    TH1F* hNstrawHits_;
    TH1F* hNcrystalHits_;

    // Tools to fill other histograms.
    DiagnosticsG4              diagnostics_;
    GeneratorSummaryHistograms genSummary_;

  };

  MinimumHits::MinimumHits(fhicl::ParameterSet const& pset):
    generatorModuleLabel_(pset.get<string>("generatorModuleLabel")),
    g4ModuleLabel_(pset.get<string>("g4ModuleLabel")),
    strawHitMakerLabel_(pset.get<string>("strawHitMakerLabel")),
    crystalHitMakerLabel_(pset.get<string>("crystalHitMakerLabel")),
    trackerStepPoints_(pset.get<string>("trackerStepPoints","tracker")),
    caloStepPoints_(pset.get<string>("caloStepPoints","calorimeter")),
    caloROStepPoints_(pset.get<string>("caloROStepPoints","calorimeterRO")),
    foilStepPoints_(pset.get<string>("foilStepPoints","stoppingtarget")),
    crvStepPoints_(pset.get<string>("CRVStepPoints","CRV")),
    vDetStepPoints_(pset.get<string>("vDetStepPoints","virtualdetector")),
    hNstrawHits_(0),
    hNcrystalHits_(0),
    diagnostics_(),
    genSummary_(){
  }

  bool MinimumHits::beginRun(art::Run& ){

    // Book histograms; must wait until beginRun because some histogram
    // limits are set using the GeometryService.
    diagnostics_.book("G4Summary");
    genSummary_.book("GeneratorSummary");

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir( "HitSummary" );

    hNstrawHits_   = tfdir.make<TH1F>( "hNstrawHits",   "Number of Straw Hits",    200, 0., 200.  );
    hNcrystalHits_ = tfdir.make<TH1F>( "hNcrystalHits", "Number of Crystal Hits",   50, 0.,  50.  );

    return true;
  }

  bool
  MinimumHits::filter(art::Event& event) {

    art::Handle<StatusG4> g4StatusHandle;
    event.getByLabel( g4ModuleLabel_, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;

    // Accept only events with good status from G4.
    if ( g4Status.status() > 1 ) {

      // Diagnostics for rejected events.
      diagnostics_.fillStatus( g4Status);
      return false;
    }

    // Get enough information to make the filter decision.
    art::Handle<StepPointMCCollection> trackerStepsHandle;
    event.getByLabel(g4ModuleLabel_, trackerStepPoints_,trackerStepsHandle);
    StepPointMCCollection const& trackerSteps(*trackerStepsHandle);

    art::Handle<StepPointMCCollection> caloStepsHandle;
    event.getByLabel(g4ModuleLabel_, caloStepPoints_, caloStepsHandle);
    StepPointMCCollection const& caloSteps(*caloStepsHandle);

    art::Handle<StepPointMCCollection> caloROStepsHandle;
    event.getByLabel(g4ModuleLabel_, caloROStepPoints_, caloROStepsHandle);
    StepPointMCCollection const& caloROSteps(*caloROStepsHandle);

    if ( trackerSteps.empty() && caloSteps.empty() && caloROSteps.empty() ){
      return false;
    }

    // Get the remaining data products from the event.
    art::Handle<GenParticleCollection> gensHandle;
    event.getByLabel( generatorModuleLabel_, gensHandle);
    GenParticleCollection const& gens(*gensHandle);

    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel(g4ModuleLabel_,simsHandle);
    SimParticleCollection const& sims(*simsHandle);

    art::Handle<StepPointMCCollection> foilStepsHandle;
    event.getByLabel(g4ModuleLabel_, foilStepPoints_, foilStepsHandle);
    StepPointMCCollection const& foilSteps(*foilStepsHandle);

    art::Handle<StepPointMCCollection> crvStepsHandle;
    event.getByLabel(g4ModuleLabel_, crvStepPoints_, crvStepsHandle);
    StepPointMCCollection const& crvSteps(*crvStepsHandle);

    art::Handle<StepPointMCCollection> vDetStepsHandle;
    event.getByLabel(g4ModuleLabel_, vDetStepPoints_, vDetStepsHandle);
    StepPointMCCollection const& vDetSteps(*vDetStepsHandle);

    art::Handle<PointTrajectoryCollection> trajectoriesHandle;
    event.getByLabel(g4ModuleLabel_,trajectoriesHandle);
    PointTrajectoryCollection const& trajectories(*trajectoriesHandle);

    art::Handle<StrawHitCollection> strawHitsHandle;
    event.getByLabel(strawHitMakerLabel_, strawHitsHandle);
    StrawHitCollection const& strawHits(*strawHitsHandle);

    art::Handle<CaloCrystalHitCollection> crystalHitsHandle;
    event.getByLabel(crystalHitMakerLabel_, crystalHitsHandle);
    CaloCrystalHitCollection const& crystalHits(*crystalHitsHandle);

    art::Handle<PhysicalVolumeInfoCollection> volsHandle;
    event.getRun().getByLabel(g4ModuleLabel_,volsHandle);
    PhysicalVolumeInfoCollection const& vols(*volsHandle);

    // Fill histograms for the events that pass the filter.
    diagnostics_.fill( g4Status,
                       sims,
                       trackerSteps,
                       caloSteps,
                       caloROSteps,
                       crvSteps,
                       foilSteps,
                       vDetSteps,
                       trajectories,
                       vols);

    genSummary_.fill( gens );

    hNstrawHits_  ->Fill( strawHits.size() );
    hNcrystalHits_->Fill( crystalHits.size() );

    return true;

  } // end of ::analyze.

}

using mu2e::MinimumHits;
DEFINE_ART_MODULE(MinimumHits);
