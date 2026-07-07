// requires a minimum number of steps in online straws

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/TrackerStatus.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "CLHEP/Random/RandFlat.h"

#include <string>

using namespace std;

namespace mu2e{

  class StationStepSelector : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;

    struct Config {
      fhicl::Atom<art::InputTag> stepstag { Name("StrawGasStepTag"), Comment("StrawGasStep tag")};
      fhicl::Atom<unsigned> minsteps { Name("MinSteps"), Comment("Minimum number of SGS")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit StationStepSelector(const Parameters& config);
    virtual void produce(art::Event& e);

  private:
    ProditionsHandle<TrackerStatus> _trackerStatus_h;
    art::RandomNumberGenerator::base_engine_t& _engine;
    CLHEP::RandFlat _randflat;

    art::InputTag _stepsTag;
    size_t        _minSteps;
  };

  StationStepSelector::StationStepSelector(const Parameters& config) :
    art::EDProducer{config},
    _engine(createEngine( art::ServiceHandle<SeedService>()->getSeed())),
    _randflat( _engine ),
    _stepsTag(config().stepstag()),
    _minSteps(config().minsteps())
    {
      produces<StrawGasStepCollection>();
    }

  void StationStepSelector::produce(art::Event& event)  {
    unique_ptr<StrawGasStepCollection> outsteps(new StrawGasStepCollection);
    auto const& trackerStatus = _trackerStatus_h.getPtr(event.id());

    std::vector<size_t> counts(StrawId::_nstations,0);
    auto steps = event.getValidHandle<StrawGasStepCollection>(_stepsTag);
    auto const& stepcol = *steps;
    for (auto const& step: stepcol){
      if (!trackerStatus->noSignal(step.strawId()))
        counts[step.strawId().station()]++;
    }
    std::vector<size_t> goodstations;
    for (size_t i=0;i<StrawId::_nstations;i++){
      if (counts[i] >= _minSteps)
        goodstations.push_back(i);
    }
    if (goodstations.size() > 0){
      size_t randStation = _randflat.fireInt(goodstations.size());
      for (auto const& step: stepcol){
        if (step.strawId().station() == goodstations[randStation])
          outsteps->push_back(step);
      }
    }
    event.put(move(outsteps));
  }
}

DEFINE_ART_MODULE(mu2e::StationStepSelector)
