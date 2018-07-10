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
      fhicl::Atom<double> mean{Name("mean"), Comment("Mean number of protons per microbunch") };
      fhicl::Atom<double> sigma{Name("sigma"), Comment("sigma of the lognormal distribution")};
      fhicl::Atom<double> cutoff{Name("cutoff"), Comment("The high tail of the distribution will be truncated at the cutoff.")};
    };

    typedef art::EDProducer::Table<Config> Parameters;

    explicit ProtonBunchIntensityLogNormal(const Parameters& conf);
    void produce(art::Event& evt) override;

  private:
    artURBG urbg_;
    std::lognormal_distribution<double> lognd_;
    double mean_;
    double cutoff_;

    static double solveForMu(double mean, double sigma, double cutoff);
    static double muMapping(double mean, double mu, double sigma, double cutoff);
  };

  ProtonBunchIntensityLogNormal::ProtonBunchIntensityLogNormal(const Parameters& conf)
    : urbg_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , lognd_(solveForMu(conf().mean(), conf().sigma(), conf().cutoff()),
             conf().sigma())
    , mean_(conf().mean())
    , cutoff_(conf().cutoff())
  {
    produces<mu2e::ProtonBunchIntensity>();
  }

  void ProtonBunchIntensityLogNormal::produce(art::Event& event) {
    double res=0.;

    while( (res = lognd_(urbg_)) > cutoff_) {}

    // convert to nearest ingeger and write out
    event.put(std::make_unique<ProtonBunchIntensity>(unsigned(rint(res)), mean_));
  }

  //================================================================
  double ProtonBunchIntensityLogNormal::solveForMu(double mean, double sigma, double cutoff) {

    if(mean <= 0.) {
      throw cet::exception("BADCONFIG")<<"ProtonBunchIntensityLogNormal: illegal mean = "
                                       <<mean<<" <= 0.\n";
    }
    if(sigma <= 0.) {
      throw cet::exception("BADCONFIG")<<"ProtonBunchIntensityLogNormal: illegal sigma = "
                                       <<sigma<<" <= 0.\n";
    }

    const double sigma2 = std::pow(sigma,2);
    double mu0 = log(mean) - 0.5*sigma2;

    // standard deviation of a non-truncated distribution (mu0,sigma)
    double sd0 = exp(mu0 +sigma2/2.)*sqrt(exp(sigma2) - 1.);

    // Our approach to solve for mu to get the requested mean works
    // for a cutoff reasonably far from the center of the distribution.
    if(cutoff < mean + sd0) {
      throw cet::exception("BADCONFIG")<<"ProtonBunchIntensityLogNormal: cutoff = "<<cutoff
                                       <<" is not above the bulk of the distribution:"
                                       <<" mean + one standard deviation ~= "<<mean + sd0
                                       <<"\n";
    }

    // We write the transcendental equation for mu as
    //
    //    mu = F(mean,mu,sigma,cutoff)
    //
    // F() here happens to be a contraction mapping, at least in the
    // range of parameters we care about.  The root can be found as
    // the stationary point of the mapping.

    const double precision = 1.e-9;
    double mu = mu0;
    for(int i=0; i<20; ++i) {
      const double tmp = muMapping(mean, mu, sigma, cutoff);
      const double delta = tmp-mu;
      mu = tmp;
      if(std::abs(delta) < precision) {
        return mu;
      }
    }

    throw cet::exception("DEBUG")<<"ProtonBunchIntensityLogNormal(): solver for mu did not converge.\n";
  }

  namespace {
    // CDF of the standard normal distribution
    double Phi(double x) {
      return 0.5*erfc(-x/sqrt(2.));
    }
  }

  double ProtonBunchIntensityLogNormal::muMapping(double mean, double mu, double sigma, double cutoff) {
    const double sigma2 = std::pow(sigma,2);
    const double argd = (log(cutoff) - mu)/sigma;
    const double res = log(mean) - sigma2/2. + log(Phi(argd)/Phi(argd - sigma));
    return res;
  }

}

DEFINE_ART_MODULE(mu2e::ProtonBunchIntensityLogNormal);
