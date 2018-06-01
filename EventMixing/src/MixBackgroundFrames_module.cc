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

  namespace {
    using namespace fhicl;
    struct MyTopConfig {
      Table<Mu2eProductMixer::Config> products { Name("products") };

      Atom<art::InputTag> protonBunchIntensityTag { Name("protonBunchIntensityTag"),
          Comment("InputTag of a ProtonBunchIntensity product representing beam fluctuations.")
          };

      Atom<double> meanEventsPerProton { Name("meanEventsPerProton"),
          Comment("The mean number of secondary events to mix per proton on target. "
                  "The number of protons on target for each output microbunch-event will be taken "
                  "from the protonBunchIntensity input."
                  ),
          1u
          };
    };

    // The following hack will hopefully go away after
    // https://cdcvs.fnal.gov/redmine/issues/19970
    // is resolved.
    MyTopConfig
    retrieveConfiguration(const std::string& subTableName, const fhicl::ParameterSet& pset)
    {
      std::set<std::string> ignorable_keys;

      // Ignore everything but the subtable
      const auto& allnames = pset.get_names();
      for(const auto& i: allnames) {
        if(i != subTableName) {
          ignorable_keys.insert(i);
        }
      }

      return fhicl::Table<MyTopConfig>(pset.get<fhicl::ParameterSet>(subTableName),
                                       ignorable_keys )();
    }
  }

  //----------------------------------------------------------------
  // Our "detail" class for art/Framework/Modules/MixFilter.h
  class MixBackgroundFramesDetail {
    Mu2eProductMixer spm_;
    art::InputTag pbiTag_;
    const double meanEventsPerProton_;
    artURBG urbg_;

    ProtonBunchIntensity pbi_;

    // intended to be called from the constructor, thus static
    static art::RandomNumberGenerator::base_engine_t& createArtEngine();

  public:
    MixBackgroundFramesDetail(const fhicl::ParameterSet& pset, art::MixHelper &helper);

    void startEvent(const art::Event& event);

    size_t nSecondaries();
  };

  //================================================================
  art::RandomNumberGenerator::base_engine_t&
  MixBackgroundFramesDetail::createArtEngine() {
    auto& engine = art::ServiceHandle<art::RandomNumberGenerator>()->getEngine();
    int dummy(0);
    engine.setSeed( art::ServiceHandle<SeedService>()->getSeed(), dummy );
    return engine;
  }

  //================================================================
  MixBackgroundFramesDetail::MixBackgroundFramesDetail(const fhicl::ParameterSet& pset, art::MixHelper& helper)
    : spm_{ retrieveConfiguration("mu2e", pset).products(), helper }
    , pbiTag_{ retrieveConfiguration("mu2e", pset).protonBunchIntensityTag() }
    , meanEventsPerProton_{ retrieveConfiguration("mu2e", pset).meanEventsPerProton() }
    , urbg_{ createArtEngine() }
  {}

  //================================================================
  void MixBackgroundFramesDetail::startEvent(const art::Event& event) {
    pbi_ = *event.getValidHandle<ProtonBunchIntensity>(pbiTag_);
  }

  //================================================================
  size_t MixBackgroundFramesDetail::nSecondaries() {
    double mean = meanEventsPerProton_ * pbi_.intensity();
    std::poisson_distribution<size_t> poisson(mean);
    auto res = poisson(urbg_);
    return res;
  }

  //================================================================
  // This is the module class.
  typedef art::MixFilter<MixBackgroundFramesDetail> MixBackgroundFrames;
}

DEFINE_ART_MODULE(mu2e::MixBackgroundFrames);
