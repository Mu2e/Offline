//
// Select events with a minimum number of StepPointMC's in various detectors.
// $Id: RPCFilter1_module.cc,v 1.3 2013/05/30 18:40:35 gandr Exp $
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
#include "TNtuple.h"

// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <iostream>

using namespace std;

namespace mu2e {

  class RPCFilter1 : public art::EDFilter {
  public:
    explicit RPCFilter1(fhicl::ParameterSet const& pset);
    virtual ~RPCFilter1() { }

    bool filter( art::Event& event);

    virtual bool beginRun(art::Run &run);
    virtual void endJob();

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
    std::string extMonUCITofStepPoints_;

    // Histogram pointers.
    TH1F* hNstrawHits_;
    TH1F* hNcrystalHits_;
    TH1F* hEDep_;
    TH1F* hEDepWide_;
    TH1F* hEDepStep_;
    TNtuple* ntup_;

    // Tools to fill other histograms.
    DiagnosticsG4              diagnostics_;
    GeneratorSummaryHistograms genSummary_;

    // Number of events that pass the filter.
    int nPassed_;

  };

  RPCFilter1::RPCFilter1(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
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
    extMonUCITofStepPoints_(pset.get<string>("extMonUCITofStepPoints","extinctionmonitorTof")),
    hNstrawHits_(0),
    hNcrystalHits_(0),
    hEDep_(0),
    hEDepWide_(0),
    hEDepStep_(0),
    ntup_(0),
    diagnostics_(),
    genSummary_(),
    nPassed_(0){
  }

  bool RPCFilter1::beginRun(art::Run& ){

    // Book histograms; must wait until beginRun because some histogram
    // limits are set using the GeometryService.
    diagnostics_.book("G4Summary");
    genSummary_.book("GeneratorSummary");

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir( "HitSummary" );

    hNstrawHits_   = tfdir.make<TH1F>( "hNstrawHits",   "Number of Straw Hits",    200, 1., 201.  );
    hNcrystalHits_ = tfdir.make<TH1F>( "hNcrystalHits", "Number of Crystal Hits",   50, 1.,  51.  );
    hEDep_         = tfdir.make<TH1F>( "hEDep",     "Energy deposition, Straw",   100, 0.,  10.  );
    hEDepWide_     = tfdir.make<TH1F>( "hEDepWide", "Energy deposition, Straw",   100, 0.,  1000. );
    hEDepStep_     = tfdir.make<TH1F>( "hEDepStep", "Energy deposition, Step",   100, 0.,  20.  );

    ntup_ = tfdir.make<TNtuple>( "ntup", "Event Summary", "x:y:z:r:p:pt:cz:ntrkSteps:nCaloSteps:nstraws:nxstals");

    return true;
  }

  bool
  RPCFilter1::filter(art::Event& event) {

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

    // Anticipate future use that selects on other combinations.
    if ( trackerSteps.empty() ){
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

    art::Handle<StepPointMCCollection> extMonUCITofStepsHandle;
    event.getByLabel(g4ModuleLabel_, extMonUCITofStepPoints_, extMonUCITofStepsHandle);
    StepPointMCCollection const& extMonUCITofSteps(*extMonUCITofStepsHandle);

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
                       extMonUCITofSteps,
                       trajectories,
                       vols);

    genSummary_.fill( gens );

    hNstrawHits_  ->Fill( strawHits.size() );
    hNcrystalHits_->Fill( crystalHits.size() );

    for ( size_t i=0; i< strawHits.size(); ++i ){
      StrawHit const& s = strawHits.at(i);
      hEDep_->Fill( s.energyDep()*1000. );
      hEDepWide_->Fill( s.energyDep()*1000. );
    }

    for ( size_t i=0; i< trackerSteps.size(); ++i ){
      StepPointMC const& s = trackerSteps.at(i);
      hEDepStep_->Fill( s.eDep()*1000. );
    }

    // Properties of generated particles, correlated with step point counts.
    float nt[ntup_->GetNvar()];
    nt[7]  = trackerSteps.size();
    nt[8]  = caloSteps.size() + caloROSteps.size();
    nt[9]  = strawHits.size();
    nt[10] = crystalHits.size();
    for ( GenParticleCollection::const_iterator i=gens.begin(), e=gens.end();
          i != e; ++i ){

      GenParticle const& gen(*i);

      CLHEP::Hep3Vector const&       pos(gen.position());
      CLHEP::Hep3Vector const&       p(gen.momentum().vect());

      double r(pos.perp());
      double cz = p.cosTheta();

      nt[0] = pos.x();
      nt[1] = pos.y();
      nt[2] = pos.z();
      nt[3] = r;
      nt[4] = p.mag();
      nt[5] = p.perp();
      nt[6] = cz;
      ntup_->Fill(nt);

    }

    nPassed_++;
    return true;

  } // end of ::analyze.

  void RPCFilter1::endJob() {
    mf::LogInfo("Summary") 
      << "RPCFilter1_module: Number of events passing the filter: " 
      << nPassed_
      << "\n";
  }

}

using mu2e::RPCFilter1;
DEFINE_ART_MODULE(RPCFilter1);
