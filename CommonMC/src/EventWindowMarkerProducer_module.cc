// ======================================================================
//
// EventWindowMarkerProducer_module:  Provides the start of event window marker
//   generated from the RF-0 signal synchronized to the system 40 MHz clock
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "SeedService/inc/SeedService.hh"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"

#include "DataProducts/inc/EventWindowMarker.hh"

#include "CLHEP/Random/RandFlat.h"

#include <iostream>
#include <string>

using namespace std;
namespace mu2e {

  class EventWindowMarkerProducer : public art::EDProducer {
    public:
      explicit EventWindowMarkerProducer(fhicl::ParameterSet const& pset);

    private:
      void beginJob() override;
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;

      float _systemClockSpeed;
      float _eventWindowConstantOffset;
      art::RandomNumberGenerator::base_engine_t& _engine;
      CLHEP::RandFlat _randflat;
  };

  EventWindowMarkerProducer::EventWindowMarkerProducer(fhicl::ParameterSet const& pset):
    art::EDProducer{pset},
    _systemClockSpeed(pset.get<float>("SystemClockSpeed",40.0)), // MHz
    _eventWindowConstantOffset(pset.get<float>("EventWindowConstantOffset",0.0)), // ns (is a calibration value to get data/sim agreement)
    // a positive value means the event window is delayed relative to the StepPointMC timing
    _engine(createEngine( art::ServiceHandle<SeedService>()->getSeed())),
    _randflat( _engine )
  {
    produces<EventWindowMarker>();
  }

  void EventWindowMarkerProducer::beginJob() {
  }

  void EventWindowMarkerProducer::beginRun(art::Run& run) {
  }

  void EventWindowMarkerProducer::produce(art::Event& event) {
      unique_ptr<EventWindowMarker> marker(new EventWindowMarker);
      marker.get()->_timeOffset = _eventWindowConstantOffset + _randflat.fire(1/_systemClockSpeed*1000) - 1/_systemClockSpeed*1000/2.; // ns
      event.put(std::move(marker));
  }

}

using mu2e::EventWindowMarkerProducer;
DEFINE_ART_MODULE(EventWindowMarkerProducer);
