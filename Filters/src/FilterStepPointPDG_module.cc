// Select hits from an StepPointMCCollection-s based on PDG ID and
// copy selected hits to output collections, preserving product
// instance names.
//
// Andrei Gaponenko, 2013, 2024

#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

// art includes.
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

namespace mu2e {

  //================================================================
  class FilterStepPointPDG : public art::EDFilter {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Sequence<art::InputTag> inputs {
        Name("inputs"),
        Comment("A list of StepPointMCCollection-s to process")
      };

      fhicl::Sequence<int> pdgToDrop {
        Name("pdgToDrop"),
        Comment("A list of particle types to drop while keeping everything else.  Mutually exclusive with pdgToKeep."),
        std::vector<int>()
      };

      fhicl::Sequence<int> pdgToKeep {
        Name("pdgToKeep"),
        Comment("A list of particle types to keep while dropping everything else.  Mutually exclusive with pdgToDrop."),
        std::vector<int>()
      };
    };

    using Parameters = art::EDFilter::Table<Config>;
    explicit FilterStepPointPDG(const Parameters& conf);

    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;

  private:

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
  };

  //================================================================
  FilterStepPointPDG::FilterStepPointPDG(const Parameters& conf)
    : art::EDFilter{conf}
    , inputs_{conf().inputs()}
    , numInputHits_()
    , numOutputHits_()
    , numInputEvents_()
    , numPassedEvents_()
  {
    for(const auto& i: inputs_) {
      // Coalesce same instance names from multiple input modules/processes.
      outNames_.insert(i.instance());
    }

    for(const auto& i : outNames_) {
      produces<StepPointMCCollection>(i);
    }

    for(const auto i: conf().pdgToDrop()) {
      pdgToDrop_.emplace_back(PDGCode::type(i));
    }

    for(const auto i: conf().pdgToKeep()) {
      pdgToKeep_.emplace_back(PDGCode::type(i));
    }

    if(pdgToDrop_.empty() && pdgToKeep_.empty()) {
      throw cet::exception("BADCONFIG")
        <<"FilterStepPointPDG: either pdgToDrop or pdgToKeep must be specified, neither is given.\n";
    }
    if(!pdgToDrop_.empty() && ! pdgToKeep_.empty()) {
      throw cet::exception("BADCONFIG")
        <<"FilterStepPointPDG: pdgToDrop and pdgToKeep settings are mutually exclusive,"
        " please specify only one of them.\n";
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

DEFINE_ART_MODULE(mu2e::FilterStepPointPDG)
