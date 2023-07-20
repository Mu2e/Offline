// This module "mixes" requested data products from a secondary input
// file into the current event.
//
// Andrei Gaponenko, 2018

#include "art/Framework/Principal/Event.h"
#include "art/Framework/IO/ProductMix/MixHelper.h"
#include "art/Framework/Modules/MixFilter.h"
#include "art_root_io/RootIOPolicy.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/TupleAs.h"
#include "canvas/Utilities/InputTag.h"

#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/EventMixing/inc/Mu2eProductMixer.hh"
#include "Offline/Mu2eUtilities/inc/artURBG.hh"
#include <random>

//================================================================
namespace mu2e {

  //----------------------------------------------------------------
  // Our "detail" class for art/Framework/Modules/MixFilter.h
  class ResamplingMixerDetail {
    Mu2eProductMixer spm_;
    const int debugLevel_;
    const unsigned maxEventsToSkip_;
    bool writeEventIDs_;
    art::EventIDSequence idseq_;
    art::RandomNumberGenerator::base_engine_t& engine_;
    artURBG urbg_;

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

      fhicl::Atom<unsigned> maxEventsToSkip { Name("MaxEventsToSkip"),
          Comment("Maximum number of events to skip at the beginning of the first secondary input file in sequential readMode.\n"),
          0u };

      fhicl::Atom<bool> writeEventIDs { Name("writeEventIDs"),
          Comment("Write out IDs of events on the secondary input stream."),
          true
          };
      fhicl::Atom<int> debugLevel { Name("debugLevel"),
          Comment("control the level of debug output"),
          0u
          };
     };

    struct Config {
      fhicl::Table<Mu2eConfig> mu2e { fhicl::Name("mu2e") };
    };

    using Parameters = art::MixFilterTable<Config>;
    ResamplingMixerDetail(const Parameters& pset, art::MixHelper &helper);

    size_t eventsToSkip() {
      static bool first(true);
      size_t result(0);
      if(first) {
        first = false;
        std::uniform_int_distribution<size_t> uniform(0, maxEventsToSkip_);
        result = uniform(urbg_);
      }
      if(debugLevel_ > 0) {
        std::cout << " Skipping " << result << " Secondaries " << std::endl;
      }
      return result;
    }
    size_t nSecondaries() const { return (size_t) 1; }

    void processEventIDs(const art::EventIDSequence& seq);

    void beginSubRun(const art::SubRun& sr);
    void startEvent(art::Event const& e);
    void finalizeEvent(art::Event& e);
    void endSubRun(art::SubRun& sr);

  };

  ResamplingMixerDetail::ResamplingMixerDetail(const Parameters& pars, art::MixHelper& helper)
    : spm_{ pars().mu2e().products(), helper }
    , debugLevel_{ pars().mu2e().debugLevel() }
    , maxEventsToSkip_{ pars().mu2e().maxEventsToSkip() }
    , writeEventIDs_{ pars().mu2e().writeEventIDs() }
    , engine_{helper.createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , urbg_{ engine_ }
    {
      if(writeEventIDs_) {
        helper.produces<art::EventIDSequence>();
      }
    }

  void ResamplingMixerDetail::processEventIDs(const art::EventIDSequence& seq) {

    spm_.processEventIDs(seq);

    if(writeEventIDs_) {
      idseq_ = seq;
    }
  }

  void ResamplingMixerDetail::beginSubRun(const art::SubRun& sr) {
    spm_.beginSubRun(sr);
  }

  void ResamplingMixerDetail::startEvent(art::Event const& e) {
    spm_.startEvent(e);
  }

  void ResamplingMixerDetail::finalizeEvent(art::Event& e) {
    if(writeEventIDs_) {
      auto o = std::make_unique<art::EventIDSequence>();
      o->swap(idseq_);
      e.put(std::move(o));
    }
  }

  void ResamplingMixerDetail::endSubRun(art::SubRun& sr) {
    spm_.endSubRun(sr);
  }

}

//================================================================
namespace mu2e {
  // This is the module class.
  typedef art::MixFilter<ResamplingMixerDetail,art::RootIOPolicy> ResamplingMixer;

}

DEFINE_ART_MODULE(mu2e::ResamplingMixer)
