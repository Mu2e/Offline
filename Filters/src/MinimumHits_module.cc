//
// Select events with a minimum number of StepPointMC's in various detectors.
// $Id: MinimumHits_module.cc,v 1.12 2013/05/30 18:40:35 gandr Exp $
// $Author: gandr $
// $Date: 2013/05/30 18:40:35 $
//
// Contact person Rob Kutschke.
//

// Mu2e includes.
#include "Mu2eUtilities/inc/DiagnosticsG4.hh"
#include "Mu2eUtilities/inc/GeneratorSummaryHistograms.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepFilterMode.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"

// Framework includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

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
    virtual void endJob();

  private:

    // Which data product collections do we want to look at.
    StepFilterMode mode_;

    int minnstraws_;
    double minpmom_;
    std::vector<PDGCode::type> pdgs_;

    int diag_, debug_;
    // Module label of the g4 module that made the generated particles
    std::string generatorModuleLabel_;

    // Module label of the module that made the StepPointMCCollections.
    std::string g4ModuleLabel_;

    // Instance names of the StepPointMCCollections.
    std::string trackerStepPoints_;
    std::string caloStepPoints_;
    std::string caloROStepPoints_;
    std::string foilStepPoints_;
    std::string crvStepPoints_;
    std::string vDetStepPoints_;
    std::string extMonUCITofStepPoints_;

    // Histogram pointers.
    TH1F* hNstraws_;

    // Tools to fill other histograms.
    DiagnosticsG4              g4Summary_;
    GeneratorSummaryHistograms genSummary_;

    // Number of events that pass the filter.
    int nPassed_;

    // Returns true if the event fails the filter.
    bool doesEventPass( unsigned nstraws, 
                        StepPointMCCollection const& calSteps,
                        StepPointMCCollection const& calROSteps,
                        StepPointMCCollection const& crvROSteps
                        );

    // count the number of unique straws hit in this event
    unsigned countStraws(StepPointMCCollection const& trkSteps);


  };

  MinimumHits::MinimumHits(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    mode_(StepFilterMode(pset.get<string>("mode"))),
    minnstraws_(pset.get<int>("MinNStraws")),
    minpmom_(pset.get<double>("MinPartMom")),
    diag_(pset.get<int>("diagLevel",0)),
    debug_(pset.get<int>("debugLevel",0)),
    generatorModuleLabel_(pset.get<string>("generatorModuleLabel")),
    g4ModuleLabel_(pset.get<string>("g4ModuleLabel")),
    trackerStepPoints_(pset.get<string>("trackerStepPoints","tracker")),
    caloStepPoints_(pset.get<string>("caloStepPoints","calorimeter")),
    caloROStepPoints_(pset.get<string>("caloROStepPoints","calorimeterRO")),
    foilStepPoints_(pset.get<string>("foilStepPoints","stoppingtarget")),
    crvStepPoints_(pset.get<string>("CRVStepPoints","CRV")),
    vDetStepPoints_(pset.get<string>("vDetStepPoints","virtualdetector")),
    extMonUCITofStepPoints_(pset.get<string>("extMonUCITofStepPoints","ExtMonUCITof")),
    hNstraws_(0),
    nPassed_(0){

      auto pdgs = pset.get<std::vector<int>>("PDGCodes");
      for(auto pdg : pdgs )
	pdgs_.push_back(PDGCode::type(pdg));


      mf::LogInfo("CONFIG")
	<< "MinimumHits_module will run in the mode: " << mode_  << "\n";
    }

  bool MinimumHits::beginRun(art::Run& ){

    static int nRuns(0);
    ++nRuns;

    if ( nRuns >= 2 ){
      if ( nRuns == 2 ){
	mf::LogInfo("CONFIG")
	  << "MinimumHits_module ignores any geometry changes at run boundaries.  Hope that's OK.\n";
      }
      return 1;
    }

    // Book histograms; must wait until beginRun because some histogram
    // limits are set using the GeometryService.
    if(diag_ > 0){
      g4Summary_.book("G4Summary");
      genSummary_.book("GeneratorSummary");
      art::ServiceHandle<art::TFileService> tfs;
      hNstraws_   = tfs->make<TH1F>( "hNstraws",   "Number of Straws",    200, 1., 201.  );
    }

    return true;
  }

  // Returns true if the event fails the filter
  bool MinimumHits::doesEventPass( unsigned nstraws, 
                                   StepPointMCCollection const& calSteps,
                                   StepPointMCCollection const& calROSteps,
                                   StepPointMCCollection const& crvSteps
                                   ){
    bool fail(true);

    switch(mode_) {

      case StepFilterMode::anyDetector      :
	fail=nstraws<(size_t)minnstraws_+1 &&
	  calSteps.empty()&&
	  calROSteps.empty()&&
	  crvSteps.empty();
	break;

      case StepFilterMode::trackerOnly      :
	fail=nstraws<(size_t)minnstraws_+1;
	break;

      case StepFilterMode::calorimeterOnly  :
	fail=calSteps.empty()&&
	  calROSteps.empty();
	break;

      case StepFilterMode::CRVOnly      :
	fail=crvSteps.empty();
	break;

      case StepFilterMode::trackerOrCalorimeter:
	fail=nstraws<(size_t)minnstraws_+1 &&
	  calSteps.empty()&&
	  calROSteps.empty();
	break;

      default:
	throw cet::exception("CONFIG")
	  << "MinimumHits module has been given a mode it does not know about.\n";
    }

    return !fail;
  }


  bool
    MinimumHits::filter(art::Event& event) {

      art::Handle<StatusG4> g4StatusHandle;
      event.getByLabel( g4ModuleLabel_, g4StatusHandle);
      StatusG4 const& g4Status = *g4StatusHandle;

      // Accept only events with good status from G4.
      if ( g4Status.status() > 1 ) {

	// Diagnostics for rejected events.
	g4Summary_.fillStatus( g4Status);
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

    art::Handle<StepPointMCCollection> crvStepsHandle;
    event.getByLabel(g4ModuleLabel_, crvStepPoints_, crvStepsHandle);
    StepPointMCCollection const& crvSteps(*crvStepsHandle);

    // count the straws
    unsigned nstraws = countStraws(trackerSteps);
    // Make filter decision
    bool pass = doesEventPass( nstraws, caloSteps, caloROSteps, crvSteps );

    if(debug_ > 0)
      std::cout << " Event has " << nstraws << " straws and status " << pass << std::endl;
    if ( !pass ){
      return false;
    }

    if(diag_ > 0){
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

      art::Handle<StepPointMCCollection> vDetStepsHandle;
      event.getByLabel(g4ModuleLabel_, vDetStepPoints_, vDetStepsHandle);
      StepPointMCCollection const& vDetSteps(*vDetStepsHandle);

      art::Handle<StepPointMCCollection> extMonUCITofStepsHandle;
      event.getByLabel(g4ModuleLabel_, extMonUCITofStepPoints_, extMonUCITofStepsHandle);
      StepPointMCCollection const& extMonUCITofSteps(*extMonUCITofStepsHandle);

      art::Handle<PointTrajectoryCollection> trajectoriesHandle;
      event.getByLabel(g4ModuleLabel_,trajectoriesHandle);
      PointTrajectoryCollection const& trajectories(*trajectoriesHandle);

      art::Handle<PhysicalVolumeInfoCollection> volsHandle;
      event.getRun().getByLabel(g4ModuleLabel_,volsHandle);
      PhysicalVolumeInfoCollection const& vols(*volsHandle);

      // Fill histograms for the events that pass the filter.
      g4Summary_.fill( g4Status,
	  sims,
	  trackerSteps,
	  caloSteps,
	  caloROSteps,
	  crvSteps,
	  foilSteps,
	  vDetSteps,
	  extMonUCITofSteps,
	  trajectories,
	  vols);

      genSummary_.fill( gens );

      hNstraws_  ->Fill( nstraws );

    }
    nPassed_++;
    return true;

  } 

  void MinimumHits::endJob() {
    mf::LogInfo("Summary")
      << "MinimumHits_module: Number of events passing the filter: "
      << nPassed_
      << "\n";
  }

  unsigned MinimumHits::countStraws(StepPointMCCollection const& trkSteps) {
    std::set<StrawId> ids;
    for (auto const& step : trkSteps) {
      auto ifnd = std::find(pdgs_.begin(),pdgs_.end(),step.simParticle()->pdgId());
      if(ifnd != pdgs_.end() && step.momentum().mag() > minpmom_)
	ids.insert(step.strawId());
    }
    return ids.size();
  }


}

using mu2e::MinimumHits;
DEFINE_ART_MODULE(MinimumHits);
