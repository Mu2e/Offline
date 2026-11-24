// FamilySplitCounter: in a cosmic stage-2 job made with Geant4 secondary splitting, counts for each SubRun the
// events, those the muon-decay rule keeps, and the sum of 1/factor over their forced families, with the rules of
// FamilySplitSource (FamilySplitRules), and writes them as a FamilySplitSummary in the SubRun. Placed in the
// stage-2 path before an event selection, it lets FamilySplitSource (countsLabel) give the stripped events the
// weight they would have without the selection. Its rules must be those of the FamilySplitSource job that reads
// the output (Production/JobConfig/cosmic/prolog.fcl: CosmicSplitCounter, from CosmicSplit.splitter).
#include <memory>
#include <string>
#include <vector>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/TableFragment.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/MCDataProducts/inc/EventWeight.hh"
#include "Offline/MCDataProducts/inc/FamilySplitSummary.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/Sources/inc/FamilySplitRules.hh"

namespace mu2e {

  class FamilySplitCounter : public art::EDProducer {
  public:
    struct Config {
      fhicl::Atom<art::InputTag> simParticles{fhicl::Name("simParticles"),
          fhicl::Comment("the SimParticles FamilySplitSource reads"), art::InputTag("compressDetStepMCs")};
      fhicl::TableFragment<FamilySplitRules::Config> rules;
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit FamilySplitCounter(Parameters const& conf);
    void produce(art::Event& event) override;
    void beginSubRun(art::SubRun& sr) override;
    void endSubRun(art::SubRun& sr) override;
    void endJob() override;

  private:
    art::InputTag simsTag_;
    FamilySplitRules rules_;
    art::InputTag copierTag_;
    std::vector<FamilySplitRules::Forced> forced_;
    unsigned long nEvents_ = 0, nKept_ = 0;
    double sumWeights_ = 0.;
    unsigned long nEventsJob_ = 0, nKeptJob_ = 0;
    double sumWeightsJob_ = 0.;
  };

  FamilySplitCounter::FamilySplitCounter(Parameters const& conf)
    : art::EDProducer{conf}
    , simsTag_(conf().simParticles())
    , rules_(conf().rules()) {
    consumes<SimParticleCollection>(simsTag_);
    if (rules_.copyFactor() > 0) {
      copierTag_ = art::InputTag(rules_.copierLabel());
      consumes<EventWeight>(copierTag_);
    }
    produces<FamilySplitSummary, art::InSubRun>();
  }

  void FamilySplitCounter::beginSubRun(art::SubRun&) {
    nEvents_ = 0;
    nKept_ = 0;
    sumWeights_ = 0.;
  }

  void FamilySplitCounter::produce(art::Event& event) {
    ++nEvents_;
    auto const& simsH = event.getValidHandle<SimParticleCollection>(simsTag_);
    if (!rules_.keepByDecay(*simsH)) { return; }
    ++nKept_;
    unsigned copyFactor = 0;
    if (rules_.copyFactor() > 0) {
      copyFactor = rules_.copyFactorFromWeight(event.getValidHandle<EventWeight>(copierTag_)->weight(),
                                               copierTag_.encode() + " in event " + std::to_string(event.event()));
    }
    rules_.findForced(*simsH, simsH.id(), copyFactor, forced_);
    for (auto const& x : forced_) { sumWeights_ += 1. / x.factor; }
  }

  void FamilySplitCounter::endSubRun(art::SubRun& sr) {
    mf::LogInfo("Summary") << "FamilySplitCounter: " << sr.id() << ": " << nEvents_ << " events, " << nKept_
                           << " kept, sum of family weights " << sumWeights_ << "\n";
    sr.put(std::make_unique<FamilySplitSummary>(nEvents_, nKept_, sumWeights_), art::fullSubRun());
    nEventsJob_ += nEvents_;
    nKeptJob_ += nKept_;
    sumWeightsJob_ += sumWeights_;
  }

  void FamilySplitCounter::endJob() {
    mf::LogInfo("Summary") << "FamilySplitCounter: " << nEventsJob_ << " events, " << nKeptJob_
                           << " kept, sum of family weights " << sumWeightsJob_ << ", decay rule dropped "
                           << rules_.nDecayDropped() << "\n";
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FamilySplitCounter)
