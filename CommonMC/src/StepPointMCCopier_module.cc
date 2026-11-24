//
// StepPointMCCopier: copy a StepPointMC collection, writing each selected step `copies` (N) times in a row
// and every other step once. Used in cosmic stage 2 (Production/JobConfig/cosmic/S2ResamplerSplit.fcl):
// Mu2eG4 then makes N independent Geant4 primaries of each stage-1 e+-/gamma above a momentum cut, and
// FamilySplitSource (copiedPrimaries) gives each one its own event with weight 1/N. The selection must be
// the same in both: the same PDG ids and momentum cut; a step is copied if |p| > minMomentum.
//
// An event with more than `maxCopied` selected steps (a shower) is not copied: every step is written once.
// Its particles then count at weight 1, so the sample stays unbiased; it only gains nothing there. This
// bounds the number of families FamilySplitSource makes from one event (at most maxCopied * N).
//
// The output has the instance name of the input, so only the label changes in Mu2eG4 inputs.primaryTag.
// The copies are identical, including their art::Ptr to the stage-1 SimParticle.
//
// The module also writes an EventWeight (no instance name) holding 1/N for an event that was copied (or
// had nothing to copy) and 1 for one that was not. FamilySplitSource reads it to tell the two apart, since
// the stage-2 compression drops the copies that leave nothing in the detector. It is for that job only: it
// is not a weight to apply in an analysis, and FamilySplitSource does not pass it on.
//

#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <string>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/MCDataProducts/inc/EventWeight.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

namespace mu2e {

  class StepPointMCCopier : public art::EDProducer {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> input{Name("input"), Comment("StepPointMC collection")};
      fhicl::Sequence<int> pdgIds{Name("pdgIds"), Comment("PDG ids (of the step's SimParticle) to copy")};
      fhicl::Atom<double> minMomentum{Name("minMomentum"), Comment("[MeV/c] steps with |p| above this are copied")};
      fhicl::Atom<unsigned> copies{Name("copies"), Comment("N >= 2: how many times a selected step is written")};
      fhicl::Atom<unsigned> maxCopied{Name("maxCopied"),
          Comment("an event with more selected steps than this is not copied; 0: no limit")};
    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit StepPointMCCopier(Parameters const& conf);
    void produce(art::Event& event) override;
    void endJob() override;

  private:
    art::InputTag input_;
    std::set<int> pdgIds_;
    double minMomentum_;
    unsigned copies_;
    unsigned maxCopied_;
    unsigned long nSteps_ = 0, nCopied_ = 0, nEvents_ = 0, nEventsWithCopies_ = 0;
    unsigned long nCapped_ = 0, nCappedSteps_ = 0;
  };

  StepPointMCCopier::StepPointMCCopier(Parameters const& conf)
    : art::EDProducer{conf}
    , input_(conf().input())
    , minMomentum_(conf().minMomentum())
    , copies_(conf().copies())
    , maxCopied_(conf().maxCopied()) {
    for (int id : conf().pdgIds()) { pdgIds_.insert(id); }
    if (copies_ < 2 || pdgIds_.empty() || minMomentum_ < 0.) {
      throw cet::exception("CONFIG") << "StepPointMCCopier: need copies >= 2, pdgIds and minMomentum >= 0\n";
    }
    consumes<StepPointMCCollection>(input_);
    produces<StepPointMCCollection>(input_.instance());
    produces<EventWeight>();
  }

  void StepPointMCCopier::produce(art::Event& event) {
    auto const& steps = *event.getValidHandle<StepPointMCCollection>(input_);
    std::vector<char> selected(steps.size(), 0);
    unsigned long nSel = 0;
    for (size_t i = 0; i < steps.size(); ++i) {
      auto const& s = steps[i];
      if (pdgIds_.count(s.simParticle()->pdgId()) == 0) { continue; }
      const double p = s.momentum().mag();
      // FamilySplitSource applies the same cut to the SimParticle that Geant4 makes of this step, whose
      // momentum agrees to rounding; refuse a step that could fall on the other side of it there
      if (std::abs(p - minMomentum_) < 1.e-6) {
        throw cet::exception("StepPointMCCopier") << "step momentum " << p << " MeV/c within 1e-6 of minMomentum\n";
      }
      if (p > minMomentum_) { selected[i] = 1; ++nSel; }
    }
    const bool capped = maxCopied_ > 0 && nSel > maxCopied_;
    const unsigned n = capped ? 1 : copies_;
    auto out = std::make_unique<StepPointMCCollection>();
    out->reserve(steps.size() + nSel * (n - 1));
    for (size_t i = 0; i < steps.size(); ++i) {
      const unsigned k = selected[i] ? n : 1;
      for (unsigned c = 0; c < k; ++c) { out->push_back(steps[i]); }
    }
    nSteps_ += steps.size();
    ++nEvents_;
    if (capped) {
      ++nCapped_;
      nCappedSteps_ += nSel;
    } else if (nSel > 0) {
      ++nEventsWithCopies_;
      nCopied_ += nSel;
    }
    event.put(std::move(out), input_.instance());
    event.put(std::make_unique<EventWeight>(1. / static_cast<double>(n)));
  }

  void StepPointMCCopier::endJob() {
    std::cout << "StepPointMCCopier " << moduleDescription().moduleLabel() << ": " << nEvents_ << " events, "
              << nSteps_ << " steps, " << nCopied_ << " written " << copies_ << " times, in " << nEventsWithCopies_
              << " events; " << nCapped_ << " events with more than " << maxCopied_ << " selected steps not copied ("
              << nCappedSteps_ << " steps)" << std::endl;
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StepPointMCCopier)
