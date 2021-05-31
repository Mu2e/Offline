// Copy step points passing time cuts to a new collection.
// Note this is a producer module, not a filter, as it is
// intended to be followed by FilterG4Out, which reads other
// StepPoint collections as well.  Having this as a filter
// could abort the trigger path prematurely.
//
// Andrei Gaponenko, 2014

#include <string>
#include <vector>
#include <algorithm>
#include <limits>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

namespace mu2e {

  //================================================================
  class SelectStepPointsByTime : public art::EDProducer {
    art::InputTag input_;
    std::string outInstanceName_;
    double cutTimeMin_;
    double cutTimeMax_;

    // statistics
    unsigned numInputHits_;
    unsigned numOutputHits_;

  public:
    explicit SelectStepPointsByTime(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event) override;
    virtual void endJob() override;
  };

  //================================================================
  SelectStepPointsByTime::SelectStepPointsByTime(const fhicl::ParameterSet& pset)
    : art::EDProducer{pset}
    , input_(pset.get<std::string>("input"))
    , outInstanceName_(pset.get<std::string>("outInstanceName"))
    , cutTimeMin_(pset.get<double>("cutTimeMin",std::numeric_limits<double>::min()))
    , cutTimeMax_(pset.get<double>("cutTimeMax",std::numeric_limits<double>::max()))
    , numInputHits_()
    , numOutputHits_()
  {
    produces<StepPointMCCollection>(outInstanceName_);
    mf::LogInfo("Info")<<"SelectStepPointsByTime: cuts for "
                       <<input_
                       <<" are: min = "<<cutTimeMin_
                       <<", max = "<<cutTimeMax_
                       <<"\n";
  }

  //================================================================
  void SelectStepPointsByTime::produce(art::Event& event) {
    std::unique_ptr<StepPointMCCollection> out(new StepPointMCCollection());

    auto ih = event.getValidHandle<StepPointMCCollection>(input_);

    for(const auto& hit : *ih) {
      if( (cutTimeMin_ <= hit.time()) && (hit.time() < cutTimeMax_)) {
        out->emplace_back(hit);
      }
    }

    numInputHits_ += ih->size();
    numOutputHits_ += out->size();

    event.put(std::move(out), outInstanceName_);
  }

  //================================================================
  void SelectStepPointsByTime::endJob() {
    mf::LogInfo("Summary")<<"SelectStepPointsByTime: passed "
                          <<numOutputHits_ <<" / "<<numInputHits_
                          <<" StepPointMCs\n";
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SelectStepPointsByTime);
