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
#include <numbers>
//#include "gsl/gsl_sf_erf.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"

#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/Mu2eUtilities/inc/artURBG.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include <boost/math/special_functions/erf.hpp>

namespace mu2e {

  class ProtonBunchIntensityLogNormal : public art::EDProducer {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<int> debug{Name("debugLevel"), Comment("debug level"),0};
      fhicl::Atom<double> extendedMean{Name("extendedMean"), Comment("Mean number of protons per microbunch for distribution without cuts") };
      fhicl::OptionalAtom<double> sigma{Name("sigma"), Comment("sigma of the lognormal distribution")};
      fhicl::OptionalAtom<double> SDF{Name("SDF"), Comment("Spill Duty Factor. sigma =sqrt(-log(SDF))")};
      fhicl::Atom<double> cutMin{Name("cutMin"), Comment("The min number of protons to generate."), 0.};
      fhicl::Atom<double> cutMax{Name("cutMax"), Comment("The high tail of the distribution will be truncated at cutMax.")};
      fhicl::Atom<art::InputTag> PP{ Name("PrimaryParticle"), Comment("PrimaryParticle")};
   };

    typedef art::EDProducer::Table<Config> Parameters;

    explicit ProtonBunchIntensityLogNormal(const Parameters& conf);
    void produce(art::Event& evt) override;
    void beginSubRun(art::SubRun & subrun ) override;

  private:
    art::InputTag pp_;
    artURBG urbg_;
    std::lognormal_distribution<double> lognd_;
    std::uniform_real_distribution<double> unitflatd_;
    int debug_;
    double mean_;
    double cutMin_;
    double cutMax_;
    static double solveForMu(double mean, double sigma);
    double xLogNormalCDFInverse(double x);

  };

  ProtonBunchIntensityLogNormal::ProtonBunchIntensityLogNormal(const Parameters& conf)
    : art::EDProducer{conf}
    , pp_(conf().PP())
    , urbg_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , unitflatd_(0.0,1.0)
    , debug_(conf().debug())
    , mean_(conf().extendedMean())
    , cutMin_(conf().cutMin())
    , cutMax_(conf().cutMax())
  {
    consumes<PrimaryParticle>(pp_);
    produces<mu2e::ProtonBunchIntensity>();
    produces<mu2e::ProtonBunchIntensity,art::InSubRun>("MeanIntensity");
    double sigma(0.0), SDF(0.0);
    if(conf().sigma(sigma) && conf().SDF(SDF))
      throw cet::exception("BADCONFIG")<<"Only one of 'sigma' or 'SDF' can be specified" << std::endl;
    if((!conf().sigma(sigma)) && (!conf().SDF(SDF)))
      throw cet::exception("BADCONFIG")<<"One of 'sigma' or 'SDF' must be specified" << std::endl;
    if(conf().SDF(SDF)) sigma = sqrt(-log(SDF));
    lognd_ = std::lognormal_distribution<double>(solveForMu(conf().extendedMean(),sigma),sigma);

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
    if(debug_ > 0)
      std::cout << "LogNormal PBI with mean " << lognd_.m() << " sigma " << lognd_.s()
      << " SDF " << exp(-lognd_.s()*lognd_.s())
      << " min " << lognd_.min() << " max " << lognd_.max() << std::endl;
  }

  void ProtonBunchIntensityLogNormal::produce(art::Event& event) {
    // determine the event type
    auto pph = event.getValidHandle<PrimaryParticle>(pp_);
    auto const& pp = *pph;
    double res=cutMax_; // force to enter the loop
    if(ProcessCode::isFromProtonBeam(pp.primaryProcess())){
      // these primaries are biased in that their production is proporitional to the PBI. In this case, sample X * lognormal(X)
      while( (res < cutMin_) || (cutMax_ <= res) ) {
        res = xLogNormalCDFInverse(unitflatd_(urbg_));
      }
    } else {
      // primaries unrelated to beam protons scale as lognormal(X)
      while( (res < cutMin_) || (cutMax_ <= res) ) {
        res = lognd_(urbg_);
      }
    }

    // convert to nearest ingeger and write out
    event.put(std::make_unique<ProtonBunchIntensity>(static_cast<unsigned long long>(llrint(res))));
  }

  //================================================================
  void ProtonBunchIntensityLogNormal::beginSubRun(art::SubRun & subrun ) {
    subrun.put(std::make_unique<ProtonBunchIntensity>(static_cast<unsigned long long>(llrint(mean_))),"MeanIntensity", art::fullSubRun());
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

  double ProtonBunchIntensityLogNormal::xLogNormalCDFInverse(double x) {
    const double mu = lognd_.m();
    const double sigma = lognd_.s();
    return exp( mu + sigma*sigma -sigma*std::numbers::sqrt2_v<double>*boost::math::erfc_inv(2.0*x));
  };

}

DEFINE_ART_MODULE(mu2e::ProtonBunchIntensityLogNormal)
