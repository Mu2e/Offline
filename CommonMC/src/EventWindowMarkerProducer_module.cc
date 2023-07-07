// ======================================================================
//
// EventWindowMarkerProducer_module:  Provides the start of event window marker
//   generated from the RF-0 signal synchronized to the system 40 MHz clock
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "Offline/SeedService/inc/SeedService.hh"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/types/OptionalAtom.h"

#include "Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/DAQConditions/inc/EventTiming.hh"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

#include <iostream>
#include <string>

using namespace std;
namespace mu2e {

  class EventWindowMarkerProducer : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<unsigned> spillType { Name( "SpillType"), Comment("Whether simulating on or off spill")};
        fhicl::Atom<bool> recoFromMCTruth{ Name("RecoFromMCTruth"), Comment("Create reco PBT from MC truth"),false};
        fhicl::Atom<float> recoFromMCTruthErr{ Name("RecoFromMCTruthErr"), Comment("If creating reco PBT from MC truth, smear by this amount in ns"),0};
        fhicl::OptionalAtom<float> initialPhaseShift{ Name("InitialPhaseShift"), Comment("Initial phase shift between beam and detector clocks in fraction of one clock tick. By default is randomized each run.")};
      };

      using Parameters = art::EDProducer::Table<Config>;
      explicit EventWindowMarkerProducer(Parameters const& config);

    private:
      void beginJob() override;
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;

      EventWindowMarker::SpillType _spillType;
      bool _recoFromMCTruth;
      float _recoFromMCTruthErr;
      bool _fixInitialPhase;
      float _initialPhaseShift;
      art::RandomNumberGenerator::base_engine_t& _engine;
      CLHEP::RandGaussQ _randgauss;
      CLHEP::RandFlat _randflat;


      ProditionsHandle<EventTiming> _eventTiming_h;
  };

  EventWindowMarkerProducer::EventWindowMarkerProducer(Parameters const& config):
    art::EDProducer{config},
    _spillType(static_cast<EventWindowMarker::SpillType>( config().spillType())),
    _recoFromMCTruth( config().recoFromMCTruth()),
    _recoFromMCTruthErr( config().recoFromMCTruthErr()),
    _fixInitialPhase( false),
    _initialPhaseShift( 0),
    _engine(createEngine( art::ServiceHandle<SeedService>()->getSeed())),
    _randgauss( _engine ),
    _randflat( _engine )
    {
      _fixInitialPhase = config().initialPhaseShift(_initialPhaseShift);
      produces<EventWindowMarker>();
      produces<ProtonBunchTimeMC>();
      produces<ProtonBunchTime>();
    }

  void EventWindowMarkerProducer::beginJob() {
    // Proton bunch is 0 to 25 ns before event window marker
    if (!_fixInitialPhase){
      _initialPhaseShift = -1*_randflat.fire(1);
    }
  }

  void EventWindowMarkerProducer::beginRun(art::Run& run) {
  }

  void EventWindowMarkerProducer::produce(art::Event& event) {
    unique_ptr<EventWindowMarker> marker(new EventWindowMarker);
    unique_ptr<ProtonBunchTimeMC> pbtmc(new ProtonBunchTimeMC);
    unique_ptr<ProtonBunchTime> pbt(new ProtonBunchTime);
    marker->_spillType = _spillType;

    EventTiming const& eventTiming = _eventTiming_h.get(event.id());

    double period = (1/eventTiming.systemClockSpeed())*1000; // in ns

    if (_spillType == EventWindowMarker::SpillType::offspill){
      pbtmc->pbtime_ = 0;
      pbt->pbtime_ = 0;
      pbt->pbterr_ = 0;
      marker->_eventLength = eventTiming.offSpillLength()*period;
    }else{
      // calculate which bin we are in from event number
      int bin = event.id().event() % eventTiming.onSpillBins();
      // determine phase for this event
      // each microbunch the proton bunch gets shifted back 5 ns from the event window
      double shiftPer = -1*period/eventTiming.onSpillBins();
      double phaseShift = (_initialPhaseShift*period + bin*shiftPer);
      if (phaseShift < -1*period)
        phaseShift += period;
      // when phase shift would be > than one clock period by next microbunch
      // the event window is one clock tick shorter
      if (phaseShift < -1*period - shiftPer)
        marker->_eventLength = (eventTiming.onSpillMaxLength()-1)*period;
      else
        marker->_eventLength = eventTiming.onSpillMaxLength()*period;

      // we set pbtime to be the time of the proton bunch in terms of the DAQ clock
      // t=0 (which is determined by the event window marker arrival time)
      pbtmc->pbtime_ = -1*eventTiming.timeFromProtonsToDRMarker() + phaseShift;
      if (_recoFromMCTruth){
        pbt->pbtime_ = pbtmc.get()->pbtime_;
        if (_recoFromMCTruthErr > 0)
          pbt->pbtime_ += _randgauss.fire(0,_recoFromMCTruthErr);
        pbt->pbterr_ = _recoFromMCTruthErr;
      }else{
        pbt->pbtime_ = -1*eventTiming.timeFromProtonsToDRMarker() - period/2.;
        pbt->pbterr_ = period/sqrt(12);
      }
    }

    event.put(std::move(marker));
    event.put(std::move(pbtmc));
    event.put(std::move(pbt));
  }
}

using mu2e::EventWindowMarkerProducer;
DEFINE_ART_MODULE(EventWindowMarkerProducer)
