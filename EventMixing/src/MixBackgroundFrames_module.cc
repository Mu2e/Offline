// This module produces "background frames" that mimic Mu2e
// microbunch-events from single particle simualtion inputs.  The
// number of particles to mix is determined based on the input
// ProtonBunchIntensity object that models beam intensity
// fluctuations.  There is a random Poisson process that is on top of
// the beam intensity fluctuations, which represents the probability
// of a secondary from a given proton creating a hit in a collection
// to be mixed.  This Poisson is sampled by the module.
//
// Andrei Gaponenko, 2018

#include <random>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/IO/ProductMix/MixHelper.h"
#include "art/Framework/Modules/MixFilter.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art_root_io/RootIOPolicy.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/TupleAs.h"
#include "canvas/Utilities/InputTag.h"

#include "EventMixing/inc/Mu2eProductMixer.hh"
#include "Mu2eUtilities/inc/artURBG.hh"
#include "SeedService/inc/SeedService.hh"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "SimulationConditions/inc/SimBookkeeper.hh"
//================================================================
namespace mu2e {

  //----------------------------------------------------------------
  // Our "detail" class for art/Framework/Modules/MixFilter.h
  class MixBackgroundFramesDetail {
    Mu2eProductMixer spm_;
    art::InputTag pbiTag_;
    double meanEventsPerProton_;
    const int debugLevel_;
    art::RandomNumberGenerator::base_engine_t& engine_;
    artURBG urbg_;

    ProtonBunchIntensity pbi_;
    int totalBkgCount_;
    float skipFactor_;

    bool writeEventIDs_;
    art::EventIDSequence idseq_;

    std::vector<std::string> simStageEfficiencyTags_;
    std::vector<double> meanEventsPerPOTFactors_;
    double eff_;

    mu2e::ProditionsHandle<mu2e::SimBookkeeper> _simbookkeeperH;
    bool mixingMeanOverride_;
  public:

    struct Mu2eConfig {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<Mu2eProductMixer::Config> products { Name("products"),
          Comment("A table specifying products to be mixed.  For each supported data type\n"
                  "there is a mixingMap sequence that defines mapping of inputs to outputs.\n"
                  "Each entry in the top-level mixingMap sequence is a sequence of two strings:\n"
                  "    [ \"InputTag\", \"outputInstanceName\" ]\n"
                  "The output instance name colon \":\" is special: it means take instance name from the input tag.\n"
                  "For example, with this config:\n"
                  "   mixingMap: [ [ \"detectorFilter:tracker\", \"tracker\" ], [ \"detectorFilter:virtualdetector\", \":\" ] ]\n"
                  "the outputs will be named \"tracker\" and \"virtualdetector\"\n"
                  )
          };

      fhicl::Atom<art::InputTag> protonBunchIntensityTag { Name("protonBunchIntensityTag"),
          Comment("InputTag of a ProtonBunchIntensity product representing beam fluctuations.")
          };

      fhicl::OptionalAtom<double> meanEventsPerProton { Name("meanEventsPerProton"),
          Comment("The mean number of secondary events to mix per proton on target. "
                  "The number of protons on target for each output microbunch-event will be taken "
                  "from the protonBunchIntensity input."
                  )
          };
      fhicl::Atom<int> debugLevel { Name("debugLevel"),
          Comment("control the level of debug output"),
          0u
          };
      fhicl::Atom<float> skipFactor { Name("skipFactor"),
          Comment("mixer will skip a number of background events between 0 and this numberr multiplied by meanEventsPerProton and PBI intensity at the start of each secondary input file."),
          1
          };

      fhicl::Atom<bool> writeEventIDs { Name("writeEventIDs"),
          Comment("Write out IDs of events on the secondary input stream."),
          false
          };

      fhicl::Sequence<std::string> simStageEfficiencyTags{ Name("simStageEfficiencyTags"),
          Comment("Sequence of strings for all the previous simulation stage efficiencies that need to be included"),
          std::vector<std::string>()
          };

      fhicl::Sequence<double> meanEventsPerPOTFactors{ Name("meanEventsPerPOTFactors"),
          Comment("Sequence of double for extra numerical factors that goes into the mean events per POT"),
          std::vector<double>()
          };
    };

    // The ".mu2e" in FHICL parameters like
    // physics.filters.somemixer.mu2e.meanEventsPerProton clearly
    // separates experiment specific settings from those provided by
    // the art framework (like "somemixer.wrapFiles").
    struct Config {
      fhicl::Table<Mu2eConfig> mu2e { fhicl::Name("mu2e") };
    };

    using Parameters = art::MixFilterTable<Config>;
    explicit MixBackgroundFramesDetail(const Parameters& pars, art::MixHelper& helper);


    size_t nSecondaries();

    size_t eventsToSkip();

    void processEventIDs(const art::EventIDSequence& seq);

    void beginSubRun(const art::SubRun& sr);
    void startEvent(const art::Event& evt);
    void finalizeEvent(art::Event& e);
    void endSubRun(art::SubRun& sr);

  };

  //================================================================
  MixBackgroundFramesDetail::MixBackgroundFramesDetail(const Parameters& pars, art::MixHelper& helper)
    : spm_{ pars().mu2e().products(), helper }
    , pbiTag_{ pars().mu2e().protonBunchIntensityTag() }
    , debugLevel_{ pars().mu2e().debugLevel() }
    , engine_{helper.createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , urbg_{ engine_ }
    , totalBkgCount_(0)
    , skipFactor_{ pars().mu2e().skipFactor() }
    , writeEventIDs_{ pars().mu2e().writeEventIDs() }
    , simStageEfficiencyTags_{ pars().mu2e().simStageEfficiencyTags() }
    , meanEventsPerPOTFactors_{ pars().mu2e().meanEventsPerPOTFactors() }
    , mixingMeanOverride_(false)
  {
    if(writeEventIDs_) {
      helper.produces<art::EventIDSequence>();
    }
    if (pars().mu2e().meanEventsPerProton(meanEventsPerProton_)) {
      mixingMeanOverride_ = true;
      if (!simStageEfficiencyTags_.empty()) {
        throw cet::exception("MixBackgroundFrames") << "You have specified a number of meanEventsPerProton *and* provided a sequence of simStageEfficiencyTags. Please supply on one or the other." << std::endl;
      }
    }
  }

  //================================================================
  void MixBackgroundFramesDetail::beginSubRun(const art::SubRun& sr) {
    spm_.beginSubRun(sr);
  }

  //================================================================
  void MixBackgroundFramesDetail::endSubRun(art::SubRun& sr) {
    spm_.endSubRun(sr);
  }

  //================================================================
  void MixBackgroundFramesDetail::startEvent(const art::Event& event) {
    pbi_ = *event.getValidHandle<ProtonBunchIntensity>(pbiTag_);
    if(debugLevel_ > 0)std::cout << " Starting event mixing, Intensity = " << pbi_.intensity() << std::endl;

    eff_ = 1;

    if (!mixingMeanOverride_) {
      SimBookkeeper const& simbookkeeper = _simbookkeeperH.get(event.id());
      for (const auto& i_simStageEff : simStageEfficiencyTags_) {
        double this_eff = simbookkeeper.getEff(i_simStageEff);
        eff_ *= this_eff;

        if (debugLevel_ > 1 && !mixingMeanOverride_) {
          std::cout << " Sim Stage Efficiency (" << i_simStageEff << ") = " << this_eff << std::endl;
          std::cout << " Cumulative Total Eff = " << eff_ << std::endl;
        }
      }
      for (const auto& i_extraFactor : meanEventsPerPOTFactors_) {
        eff_ *= i_extraFactor;
        if (debugLevel_ > 1 && !mixingMeanOverride_) {
          std::cout << " Extra meanEventsPerPOT Factor = " << i_extraFactor << std::endl;
          std::cout << " Cumulative Total Eff = " << eff_ << std::endl;
        }
      }
    }
  }

  //================================================================
  size_t MixBackgroundFramesDetail::nSecondaries() {
    double mean = pbi_.intensity();
    if(mixingMeanOverride_) {
      mean *= meanEventsPerProton_;
    }
    else {
      mean *= eff_;
    }
    std::poisson_distribution<size_t> poisson(mean);
    auto res = poisson(urbg_);
    if(debugLevel_ > 0)std::cout << " Mixing " << res  << " Secondaries " << std::endl;
    return res;
  }

  //================================================================
  size_t MixBackgroundFramesDetail::eventsToSkip() {
    //FIXME: Ideally, we would know the number of events in the secondary input file
    double skipFactor = skipFactor_*pbi_.intensity();
    if(mixingMeanOverride_) {
      skipFactor *= meanEventsPerProton_;
    }
    else {
      skipFactor *= eff_;
    }
    std::uniform_int_distribution<size_t> uniform(0, skipFactor);
    size_t result = uniform(urbg_);
    if(debugLevel_ > 0) {
      std::cout << " Skipping " << result << " Secondaries " << std::endl;
    }
    return result;
  }

  //================================================================
  void MixBackgroundFramesDetail::processEventIDs(art::EventIDSequence const& seq) {
    if(writeEventIDs_) {
      idseq_ = seq;
    }

    if (debugLevel_ > 4) {
      std::cout << "The following bkg events were mixed in (START)" << std::endl;
      int counter = 0;
      for (const auto& i_eid : seq) {
        std::cout << "Run: " << i_eid.run() << " SubRun: " << i_eid.subRun() << " Event: " << i_eid.event() << std::endl;
        ++counter;
      }
      totalBkgCount_ += counter;
      std::cout << "Bkg Event Count  (this microbunch) = " << counter << std::endl;
      std::cout << "Bkg Event Count  (all microbunches) = " << totalBkgCount_ << " (END)" << std::endl;
    }
  }

  //================================================================
  void MixBackgroundFramesDetail::finalizeEvent(art::Event& e) {
    if(writeEventIDs_) {
      auto o = std::make_unique<art::EventIDSequence>();
      o->swap(idseq_);
      e.put(std::move(o));
    }
  }

  //================================================================

  //================================================================
  // This is the module class.
  typedef art::MixFilter<MixBackgroundFramesDetail,art::RootIOPolicy> MixBackgroundFrames;
}

DEFINE_ART_MODULE(mu2e::MixBackgroundFrames);
