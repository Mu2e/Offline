// requires a minimum number of steps in online straws

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/TrackerStatus.hh"

#include <string>

using namespace std;

namespace mu2e{

  class StrawGasStepFilter : public art::EDFilter {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;

    struct Config {
      fhicl::Atom<art::InputTag> stepstag { Name("StrawGasStepTag"), Comment("StrawGasStep tag")};
      fhicl::Atom<unsigned> minsteps { Name("MinSteps"), Comment("Minimum number of SGS")};
    };
    using Parameters = art::EDFilter::Table<Config>;
    explicit StrawGasStepFilter(const Parameters& config);
    virtual bool filter(art::Event& e);

  private:
    art::InputTag _stepsTag;
    size_t        _minSteps;

    ProditionsHandle<TrackerStatus> _trackerStatus_h;


  };

  StrawGasStepFilter::StrawGasStepFilter(const Parameters& config) :
    art::EDFilter{config},
    _stepsTag(config().stepstag()),
    _minSteps(config().minsteps())
    {
    }

  bool StrawGasStepFilter::filter(art::Event& event)  {
    auto const& trackerStatus = _trackerStatus_h.getPtr(event.id());

    size_t count = 0;
    auto steps = event.getValidHandle<StrawGasStepCollection>(_stepsTag);
    auto const& stepcol = *steps;
    for (auto const& step : stepcol){
      if (!trackerStatus->noSignal(step.strawId()))
        count += 1;
    }
    return count >= _minSteps;
  }
}

DEFINE_ART_MODULE(mu2e::StrawGasStepFilter)
