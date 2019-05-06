// Module to simulate the microbunch-to-microbunch intensity fluctuations of the
// protons hitting the primary target.  This module does not describe the time
// structure of the proton bunches, which is assumed to be independent of the
// intensity.  The data product of this module is used by downstream simulation
// modules, in particular, all generators which start from stopped muons need to
// use the same (coherent) intensity fluctuations for the rates of signal and
// background processes to be consistent.
//
// Original author: David Brown (LBNL) 19 May 2015
// Revamp: Andrei Gaponenko, 2018

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "SeedService/inc/SeedService.hh"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "Mu2eUtilities/inc/artURBG.hh"

#include <random>

namespace mu2e {

  class ProtonBunchIntensityFlat : public art::EDProducer {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<double> mean{Name("mean"), Comment("Mean number of protons per microbunch") };
      fhicl::Atom<double> halfWidth{Name("halfWidth"), Comment("Fractional half width of the flat distribution, 0 <= halfWidth <= 1.")};
    };

    typedef art::EDProducer::Table<Config> Parameters;

    explicit ProtonBunchIntensityFlat(const Parameters& conf);
    void produce(art::Event& evt) override;
    void beginSubRun(art::SubRun & subrun ) override;

  private:
    artURBG urbg_;
    std::uniform_real_distribution<double> flat_;
    double mean_;
  };

  ProtonBunchIntensityFlat::ProtonBunchIntensityFlat(const Parameters& conf)
    : art::EDProducer{conf}
    , urbg_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , flat_(-conf().halfWidth(), +conf().halfWidth())
    , mean_(conf().mean())
  {
    produces<mu2e::ProtonBunchIntensity>();
    produces<mu2e::ProtonBunchIntensity,art::InSubRun>("MeanIntensity");

    if(conf().halfWidth() < 0.) {
      throw cet::exception("BADCONFIG")<<"ProtonBunchIntensityFlat: illegal halfWidth = "
                                       <<conf().halfWidth()<<" < 0.\n";
    }

    if(1. < conf().halfWidth()) {
      throw cet::exception("BADCONFIG")<<"ProtonBunchIntensityFlat: illegal halfWidth = "
                                       <<conf().halfWidth()<<" > 1.\n";
    }
  }

  void ProtonBunchIntensityFlat::produce(art::Event& event) {
    double res = mean_ * (1. + flat_(urbg_));

    // convert to nearest ingeger and write out
    event.put(std::make_unique<ProtonBunchIntensity>(unsigned(rint(res))));
  }
  
  void ProtonBunchIntensityFlat::beginSubRun(art::SubRun & subrun ) {
    subrun.put(std::make_unique<ProtonBunchIntensity>(unsigned(rint(mean_))),"MeanIntensity");
  }

}

DEFINE_ART_MODULE(mu2e::ProtonBunchIntensityFlat);
