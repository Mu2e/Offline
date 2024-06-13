//  Accept events in which the specified StepPointMCCollection has
//  at least the specified number of entries.
//
// Contact person Rob Kutschke

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

#include <string>

using namespace std;

namespace mu2e{

  class STMStepPointFilter : public art::EDFilter {
  public:
    explicit STMStepPointFilter(fhicl::ParameterSet const& pset);
    virtual bool filter(art::Event& e);

  private:
    art::InputTag _stepsTag;
    size_t        _minStepPoints;

  };

  STMStepPointFilter::STMStepPointFilter(fhicl::ParameterSet const& pset) :
    art::EDFilter{pset},
    _stepsTag(pset.get<string>("stepsTag")),
    _minStepPoints(pset.get<size_t>("minStepPoints",1)){
    }

  bool STMStepPointFilter::filter(art::Event& event)  {
    auto steps = event.getValidHandle<StepPointMCCollection>(_stepsTag);
    return ( steps->size() >= _minStepPoints );
  }
}

DEFINE_ART_MODULE(mu2e::STMStepPointFilter)
