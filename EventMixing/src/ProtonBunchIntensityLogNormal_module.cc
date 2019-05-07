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

#include <random>
#include "gsl/gsl_sf_erf.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "SeedService/inc/SeedService.hh"
#include "Mu2eUtilities/inc/artURBG.hh"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"

namespace mu2e {

  class ProtonBunchIntensityLogNormal : public art::EDProducer {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<double> extendedMean{Name("extendedMean"), Comment("Mean number of protons per microbunch for distribution without cuts") };
      fhicl::Atom<double> sigma{Name("sigma"), Comment("sigma of the lognormal distribution")};
      fhicl::Atom<double> cutMin{Name("cutMin"), Comment("The min number of protons to generate."), 0.};
      fhicl::Atom<double> cutMax{Name("cutMax"), Comment("The high tail of the distribution will be truncated at cutMax.")};
    };

    typedef art::EDProducer::Table<Config> Parameters;

    explicit ProtonBunchIntensityLogNormal(const Parameters& conf);
    void produce(art::Event& evt) override;
    void beginSubRun(art::SubRun & subrun ) override;

  private:
    artURBG urbg_;
    std::lognormal_distribution<double> lognd_;
    double mean_;
    double cutMin_;
    double cutMax_;

    static double solveForMu(double mean, double sigma);
  };

  ProtonBunchIntensityLogNormal::ProtonBunchIntensityLogNormal(const Parameters& conf)
    : art::EDProducer{conf}
    , urbg_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , lognd_(solveForMu(conf().extendedMean(), conf().sigma()),
             conf().sigma())
    , mean_(conf().extendedMean())
    , cutMin_(conf().cutMin())
    , cutMax_(conf().cutMax())
  {
    produces<mu2e::ProtonBunchIntensity>();
    produces<mu2e::ProtonBunchIntensity,art::InSubRun>("MeanIntensity");

    if(cutMin_ < 0.) {
      throw cet::exception("BADCONFIG")<<"ProtonBunchIntensityLogNormal: illegal cutMin = "
                                       <<cutMin_<<" < 0.\n";
    }

    if(cutMax_ <= cutMin_) {
      throw cet::exception("BADCONFIG")<<"ProtonBunchIntensityLogNormal: illegal cutMax = "
                                       <<cutMax_<<" <= cutMin = "<<cutMin_<<"\n";
    }

    // Compute generation efficiency
    const double mu = lognd_.m();
    const double sigma = lognd_.s();
    auto logNormalCDF = [mu, sigma](double x) {
      return 0.5*erfc( - (log(x) - mu)/(sigma*sqrt(2.)) );
    };

    const double bc1 =  (0. < cutMin_) ? logNormalCDF(cutMin_) : 0.;
    const double bc2 = logNormalCDF(cutMax_);
    const double acceptance = bc2 - bc1;

    if(acceptance < 1.e-6) {
      throw cet::exception("BADCONFIG")<<"ProtonBunchIntensityLogNormal: the requested settings"
                                       <<" result in a very low generation efficiency = "<<acceptance
                                       <<".  Write a better algorithm for the task."<<"\n";
    }
  }

  void ProtonBunchIntensityLogNormal::produce(art::Event& event) {
    double res=cutMax_; // force to enter the loop
    while( (res < cutMin_) || (cutMax_ <= res) ) {
        res = lognd_(urbg_);
    }

    // convert to nearest ingeger and write out
    event.put(std::make_unique<ProtonBunchIntensity>(unsigned(rint(res))));
  }

  //================================================================
  void ProtonBunchIntensityLogNormal::beginSubRun(art::SubRun & subrun ) {
    subrun.put(std::make_unique<ProtonBunchIntensity>(unsigned(rint(mean_))),"MeanIntensity");
  }

  //================================================================
  double ProtonBunchIntensityLogNormal::solveForMu(double mean, double sigma) {
    if(mean <= 0.) {
      throw cet::exception("BADCONFIG")<<"ProtonBunchIntensityLogNormal: illegal mean = "
                                       <<mean<<" <= 0.\n";
    }
    if(sigma <= 0.) {
      throw cet::exception("BADCONFIG")<<"ProtonBunchIntensityLogNormal: illegal sigma = "
                                       <<sigma<<" <= 0.\n";
    }

    const double sigma2 = std::pow(sigma,2);

    // log-normal mu that gives distribution with the requested mean (before the cuts!)
    double mu0 = log(mean) - 0.5*sigma2;

    return mu0;
  }

}

DEFINE_ART_MODULE(mu2e::ProtonBunchIntensityLogNormal);
