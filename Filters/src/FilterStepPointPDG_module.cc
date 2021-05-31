// Select hits from an StepPointMCCollection-s based on PDG ID and
// copy selected hits to output collections, preserving product
// instance names.
//
// Andrei Gaponenko, 2013

#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "DataProducts/inc/PDGCode.hh"

namespace mu2e {

  //================================================================
  class FilterStepPointPDG : public art::EDFilter {
    typedef std::vector<art::InputTag> InputTags;
    InputTags inputs_;

    // Output instance names.
    typedef std::set<std::string> OutputNames;
    OutputNames outNames_;

    // Particle type lists are typically short, so a vector works
    // better than a set.
    std::vector<PDGCode::type> pdgToDrop_;
    std::vector<PDGCode::type> pdgToKeep_;

    // statistics
    unsigned numInputHits_;
    unsigned numOutputHits_;

    unsigned numInputEvents_;
    unsigned numPassedEvents_;

  public:
    explicit FilterStepPointPDG(const fhicl::ParameterSet& pset);
    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
  };

  //================================================================
  FilterStepPointPDG::FilterStepPointPDG(const fhicl::ParameterSet& pset)
    : art::EDFilter{pset}
    , numInputHits_()
    , numOutputHits_()
    , numInputEvents_()
    , numPassedEvents_()
  {
    const auto tags(pset.get<std::vector<std::string> >("inputs"));
    for(const auto& i : tags) {
      inputs_.emplace_back(i);
      // Coalesce same instance names from multiple input modules/processes.
      outNames_.insert(inputs_.back().instance());
    }
    for(const auto& i : outNames_) {
      produces<StepPointMCCollection>(i);
    }

    const auto drop(pset.get<std::vector<int> >("pdgToDrop", std::vector<int>()));
    for(const auto i : drop) {
      pdgToDrop_.emplace_back(PDGCode::type(i));
    }

    const auto keep(pset.get<std::vector<int> >("pdgToKeep", std::vector<int>()));
    for(const auto i : keep) {
      pdgToKeep_.emplace_back(PDGCode::type(i));
    }

    if(drop.empty() && keep.empty()) {
      throw cet::exception("BADCONFIG")
        <<"FilterStepPointPDG: either pdgToDrop or pdgToKeep must be specified.\n";
    }
    if(!drop.empty() && ! keep.empty()) {
      throw cet::exception("BADCONFIG")
        <<"FilterStepPointPDG: either pdgToDrop or pdgToKeep,"
        <<" but not both, can be used at a time.\n";
    }
  }

  //================================================================
  bool FilterStepPointPDG::filter(art::Event& event) {
    bool passed = false;

    typedef std::map<std::string, std::unique_ptr<StepPointMCCollection> > OutMap;
    OutMap outHits;
    for(const auto& i : outNames_) {
      std::unique_ptr<StepPointMCCollection> p(new StepPointMCCollection());
      outHits.insert(std::move(std::make_pair(i, std::move(p))));
    }


    for(const auto& tag : inputs_) {

      auto ih = event.getValidHandle<StepPointMCCollection>(tag);
      StepPointMCCollection& out = *outHits[tag.instance()];

      for(const auto& hit : *ih) {
        const PDGCode::type pdgId = hit.simParticle()->pdgId();

        if(!pdgToDrop_.empty() &&
           (std::find(pdgToDrop_.begin(), pdgToDrop_.end(), pdgId) == pdgToDrop_.end()))
          {
            out.emplace_back(hit);
            passed |= true;
          }

        if(!pdgToKeep_.empty() &&
           (std::find(pdgToKeep_.begin(), pdgToKeep_.end(), pdgId) != pdgToKeep_.end()))
          {
            out.emplace_back(hit);
            passed |= true;
          }
      }

      numInputHits_ += ih->size();
      numOutputHits_ += out.size();
    }

    for(const auto& i : outNames_) {
      event.put(std::move(outHits[i]), i);
    }

    ++numInputEvents_;
    if(passed) {
      ++numPassedEvents_;
    }

    return passed;
  }

  //================================================================
  void FilterStepPointPDG::endJob() {
    mf::LogInfo("Summary")<<"FilterStepPointPDG: passed "
                          <<numOutputHits_ <<" / "<<numInputHits_
                          <<" StepPointMCs, "
                          <<numPassedEvents_<<" / "<<numInputEvents_
                          <<" events\n";
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FilterStepPointPDG);
