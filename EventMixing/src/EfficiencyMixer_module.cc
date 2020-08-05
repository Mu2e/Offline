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

//================================================================
namespace mu2e {

  //----------------------------------------------------------------
  // Our "detail" class for art/Framework/Modules/MixFilter.h
  class EfficiencyMixerDetail {
    Mu2eProductMixer spm_;
    const int debugLevel_;

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

      fhicl::Atom<int> debugLevel { Name("debugLevel"),
          Comment("control the level of debug output"),
          0u
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
    explicit EfficiencyMixerDetail(const Parameters& pars, art::MixHelper& helper);

    void startEvent(const art::Event& event);

    void endRun(art::Run& r);

    size_t nSecondaries();

    void processEventIDs(const art::EventIDSequence& seq);
  };

  //================================================================
  EfficiencyMixerDetail::EfficiencyMixerDetail(const Parameters& pars, art::MixHelper& helper)
    : spm_{ pars().mu2e().products(), helper }
    , debugLevel_{ pars().mu2e().debugLevel() }
  {
    //    helper.produces<mu2e::SimEfficiencyStage, art::InRun>();
  }

  //================================================================
  void EfficiencyMixerDetail::endRun(art::Run& r) {
  }

  //================================================================
  void EfficiencyMixerDetail::startEvent(const art::Event& event) {
  }

  //================================================================
  size_t EfficiencyMixerDetail::nSecondaries() {
    return 1;
  }

  //==============================================================
  void EfficiencyMixerDetail::processEventIDs(art::EventIDSequence const& seq) {
  }


  //================================================================
  // This is the module class.
  typedef art::MixFilter<EfficiencyMixerDetail,art::RootIOPolicy> EfficiencyMixer;
}

DEFINE_ART_MODULE(mu2e::EfficiencyMixer);
