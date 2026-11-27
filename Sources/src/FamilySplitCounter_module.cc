// FamilySplitCounter: in a cosmic stage-2 job made with Geant4 secondary splitting, counts for each SubRun the
// events, those the muon-decay rule keeps, and the sum of the weights of their forced families, with the rules of
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
    FamilySplitSummary subRun_, job_;
  };

  FamilySplitCounter::FamilySplitCounter(Parameters const& conf)
    : art::EDProducer{conf}
    , simsTag_(conf().simParticles())
    , rules_(conf().rules()) {
    consumes<SimParticleCollection>(simsTag_);
    produces<FamilySplitSummary, art::InSubRun>();
  }

  void FamilySplitCounter::beginSubRun(art::SubRun&) { subRun_ = FamilySplitSummary(); }

  void FamilySplitCounter::produce(art::Event& event) {
    auto const& simsH = event.getValidHandle<SimParticleCollection>(simsTag_);
    rules_.count(*simsH, simsH.id(), subRun_);
  }

  void FamilySplitCounter::endSubRun(art::SubRun& sr) {
    mf::LogInfo("Summary") << "FamilySplitCounter: " << sr.id() << ": " << subRun_.nEvents() << " events, "
                           << subRun_.nKept() << " kept, sum of family weights " << subRun_.sumWeights() << "\n";
    sr.put(std::make_unique<FamilySplitSummary>(subRun_), art::fullSubRun());
    job_.aggregate(subRun_);
  }

  void FamilySplitCounter::endJob() {
    mf::LogInfo("Summary") << "FamilySplitCounter: " << job_.nEvents() << " events, " << job_.nKept()
                           << " kept, sum of family weights " << job_.sumWeights() << ", decay rule dropped "
                           << rules_.nDecayDropped() << "\n";
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FamilySplitCounter)
