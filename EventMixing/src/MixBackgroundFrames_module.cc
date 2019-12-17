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
      Atom<int> debugLevel { Name("debugLevel"),
          Comment("control the level of debug output"),
          0u
          };
      Atom<float> skipFactor { Name("skipFactor"),
	  Comment("mixer will skip a number of background events between 0 and this numberr multiplied by meanEventsPerProton and PBI intensity at the start of each secondary input file."),
	  1
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
    const int debugLevel_;
    art::RandomNumberGenerator::base_engine_t& engine_;
    artURBG urbg_;

    ProtonBunchIntensity pbi_;
    int totalBkgCount_;
    float skipFactor_;

  public:
    MixBackgroundFramesDetail(const fhicl::ParameterSet& pset, art::MixHelper &helper);

    void startEvent(const art::Event& event);

    size_t nSecondaries();

    size_t eventsToSkip();

    void processEventIDs(art::EventIDSequence const& seq);
  };

  //================================================================
  MixBackgroundFramesDetail::MixBackgroundFramesDetail(const fhicl::ParameterSet& pset, art::MixHelper& helper)
    : spm_{ retrieveConfiguration("mu2e", pset).products(), helper }
    , pbiTag_{ retrieveConfiguration("mu2e", pset).protonBunchIntensityTag() }
    , meanEventsPerProton_{ retrieveConfiguration("mu2e", pset).meanEventsPerProton() }
    , debugLevel_{ retrieveConfiguration("mu2e", pset).debugLevel() }
    , engine_{helper.createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , urbg_{ engine_ }
    , totalBkgCount_(0)
    , skipFactor_{ retrieveConfiguration("mu2e", pset).skipFactor() }
  {}

  //================================================================
  void MixBackgroundFramesDetail::startEvent(const art::Event& event) {
    pbi_ = *event.getValidHandle<ProtonBunchIntensity>(pbiTag_);
    if(debugLevel_ > 0)std::cout << " Starting event mixing, Intensity = " << pbi_.intensity() << std::endl;
  }

  //================================================================
  size_t MixBackgroundFramesDetail::nSecondaries() {
    double mean = meanEventsPerProton_ * pbi_.intensity();
    std::poisson_distribution<size_t> poisson(mean);
    auto res = poisson(urbg_);
    if(debugLevel_ > 0)std::cout << " Mixing " << res  << " Secondaries " << std::endl;
    return res;
  }

  //================================================================
  size_t MixBackgroundFramesDetail::eventsToSkip() {
    //FIXME: Ideally, we would know the number of events in the secondary input file
    std::uniform_int_distribution<size_t> uniform(0, skipFactor_*meanEventsPerProton_*pbi_.intensity());
    size_t result = uniform(urbg_);
    if(debugLevel_ > 0) { 
      std::cout << " Skipping " << result << " Secondaries " << std::endl; 
    }
    return result;
  }

  //================================================================
  void MixBackgroundFramesDetail::processEventIDs(art::EventIDSequence const& seq) {

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
  // This is the module class.
  typedef art::MixFilter<MixBackgroundFramesDetail,art::RootIOPolicy> MixBackgroundFrames;
}

DEFINE_ART_MODULE(mu2e::MixBackgroundFrames);
