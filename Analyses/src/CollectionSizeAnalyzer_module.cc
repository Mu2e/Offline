// Andrei Gaponenko, 2015

#include <string>
#include <vector>

#include "TDirectory.h"
#include "TProfile.h"
#include "TH1.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "MCDataProducts/inc/StepPointMCCollection.hh"

namespace mu2e {

  //================================================================
  class CollectionSizeAnalyzer : public art::EDAnalyzer {
    TProfile *hStepPointSize_;
    TH1D *hStepPointSum_;
    bool useModuleLabel_;
    bool useInstanceName_;
    bool useProcessName_;
  public:
    explicit CollectionSizeAnalyzer(const fhicl::ParameterSet& pset);
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  CollectionSizeAnalyzer::CollectionSizeAnalyzer(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , hStepPointSize_(0)
    , useModuleLabel_(pset.get<bool>("useModuleLabel"))
    , useInstanceName_(pset.get<bool>("useInstanceName"))
    , useProcessName_(pset.get<bool>("useProcessName"))
  {
    art::ServiceHandle<art::TFileService> tfs;
    hStepPointSize_ = tfs->make<TProfile>("stepPointsSize", "Average collection size", 1, 0., 1.);
    hStepPointSum_ = tfs->make<TH1D>("stepPointsSum", "Sum of step point collection entries", 1, 0., 1.);
  }

  //================================================================
  void CollectionSizeAnalyzer::analyze(const art::Event& event) {
    std::vector<art::Handle<StepPointMCCollection> > stepHandles;
    event.getManyByType(stepHandles);
    for(const auto& c: stepHandles) {
      std::ostringstream collName;

      if(useModuleLabel_) {
        collName<<c.provenance()->moduleLabel();
      }
      if(useInstanceName_) {
        if(!collName.str().empty()) {
          collName<<":";
        }
        collName<<c.provenance()->productInstanceName();
      }
      if(useProcessName_) {
        if(!collName.str().empty()) {
          collName<<":";
        }
        collName<<c.provenance()->processName();
      }

      hStepPointSize_->Fill(collName.str().c_str(), double(c->size()));
      hStepPointSum_->Fill(collName.str().c_str(), double(c->size()));
    }

  } // analyze(event)

    //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CollectionSizeAnalyzer);
