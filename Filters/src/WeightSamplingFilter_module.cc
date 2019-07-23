
// ======================================================================
//
// WeightSamplingFilterProducer_module:  Allows filtering based on eventweights
//   in order to create new subset of events with correct distribution
//
// ======================================================================

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "SeedService/inc/SeedService.hh"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"

#include "CLHEP/Random/RandFlat.h"

// Mu2e includes.
#include "MCDataProducts/inc/EventWeight.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"



#include <iostream>
#include <string>

using namespace std;
namespace mu2e {

  class WeightSamplingFilter : public art::EDFilter {
    public:
      explicit WeightSamplingFilter(fhicl::ParameterSet const& pset);

    private:
      bool beginRun(art::Run& run) override;
      bool endRun(art::Run& run) override;
      bool filter(art::Event& event) override;
      
      art::InputTag _evtWtModule;
      art::InputTag _genParticleModule;
      art::RandomNumberGenerator::base_engine_t& _engine;
      CLHEP::RandFlat _randflat;
      double _weightScalingFactor; // allows scaling weight to maximize efficiency, want max rescaled weight=>1
      int _debug;
      unsigned      _nevt, _npass;
  };

  WeightSamplingFilter::WeightSamplingFilter(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    _evtWtModule(pset.get<art::InputTag>("EventWeightModule")),
    _genParticleModule(pset.get<std::string>("genParticleModule","compressDigiMCs")),
    _engine(createEngine( art::ServiceHandle<SeedService>()->getSeed())),
    _randflat( _engine ),
    _weightScalingFactor(pset.get<double>("WeightScalingFactor",1)),
    _debug        (pset.get<int>          ("debugLevel",0)),
    _nevt(0), _npass(0){}

  bool WeightSamplingFilter::beginRun(art::Run& run) {
    return true;
  }
  bool WeightSamplingFilter::endRun(art::Run& run) {
    return true;
  }

  bool WeightSamplingFilter::filter(art::Event& event) {
    ++_nevt;

    double evtwt = event.getValidHandle<EventWeight>( _evtWtModule )->weight() * _weightScalingFactor;
    if (_debug > 0 && evtwt > 1){
      std::cout << moduleDescription().moduleLabel() << " weight scaling too high, probability is " << evtwt << " " << _weightScalingFactor << std::endl;
    }
    double temp = _randflat.fire();
    if (temp > evtwt)
      return false;

    ++_npass;
    return true;
  }
}

using mu2e::WeightSamplingFilter;
DEFINE_ART_MODULE(WeightSamplingFilter);
