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

#include "MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "RecoDataProducts/inc/ProtonBunchTime.hh"
#include "DataProducts/inc/EventWindowMarker.hh"

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
        fhicl::Atom<float> systemClockSpeed{ Name("SystemClockSpeed"), Comment("Detector clock speed in MHz"), 40.0};
        fhicl::Atom<float> constantOffset{ Name("ConstantTimeOffset"), Comment("Constant time offset correction in ns"), 0.0};
        fhicl::Atom<unsigned> spillType { Name( "SpillType"), Comment("Whether simulating on or off spill"), EventWindowMarker::SpillType::onspill};
        fhicl::Atom<unsigned> offSpillLength{ Name("OffSpillEventLength"), Comment("Length of off spill events in clock counts" ),4000};
        fhicl::Atom<int> onSpillBins{ Name("OnSpillBins"), Comment("Number of microbunches before proton bunch phase repeats"), 5};
        fhicl::Atom<unsigned> onSpillMaxLength{ Name("OnSpillEventMaxLength"), Comment("Length of longest on spill events in clock counts"),68};
        fhicl::Atom<bool> recoFromMCTruth{ Name("RecoFromMCTruth"), Comment("Create reco PBT from MC truth"),false};
        fhicl::Atom<float> recoFromMCTruthErr{ Name("RecoFromMCTruthErr"), Comment("If creating reco PBT from MC truth, smear by this amount in ns"),0};
      };

      using Parameters = art::EDProducer::Table<Config>;
      explicit EventWindowMarkerProducer(Parameters const& config);

    private:
      void beginJob() override;
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;

      float _systemClockSpeed;
      float _constantOffset;
      EventWindowMarker::SpillType _spillType;
      uint16_t _offSpillLength;
      int _onSpillBins;
      uint16_t _onSpillMaxLength;
      bool _recoFromMCTruth;
      float _recoFromMCTruthErr;
      art::RandomNumberGenerator::base_engine_t& _engine;
      CLHEP::RandGaussQ _randgauss;
      CLHEP::RandFlat _randflat;

      double _initialPhaseShift;
  };

  EventWindowMarkerProducer::EventWindowMarkerProducer(Parameters const& config):
    art::EDProducer{config},
    _systemClockSpeed( config().systemClockSpeed()), // MHz
    _constantOffset( config().constantOffset()), // ns (is a calibration value to get data/sim agreement)
    _spillType(static_cast<EventWindowMarker::SpillType>( config().spillType())),
    _offSpillLength( config().offSpillLength()),
    _onSpillBins( config().onSpillBins()),
    _onSpillMaxLength( config().onSpillMaxLength()),
    _recoFromMCTruth( config().recoFromMCTruth()),
    _recoFromMCTruthErr( config().recoFromMCTruthErr()),
    _engine(createEngine( art::ServiceHandle<SeedService>()->getSeed())),
    _randgauss( _engine ),
    _randflat( _engine )
  {
    produces<EventWindowMarker>();
    produces<ProtonBunchTimeMC>();
    produces<ProtonBunchTime>();
  }

  void EventWindowMarkerProducer::beginJob() {
    // Proton bunch is 0 to 25 ns before event window marker
    _initialPhaseShift = -1*_randflat.fire(1/_systemClockSpeed*1000);
  }

  void EventWindowMarkerProducer::beginRun(art::Run& run) {
  }

  void EventWindowMarkerProducer::produce(art::Event& event) {
      unique_ptr<EventWindowMarker> marker(new EventWindowMarker);
      unique_ptr<ProtonBunchTimeMC> pbtmc(new ProtonBunchTimeMC);
      unique_ptr<ProtonBunchTime> pbt(new ProtonBunchTime);
      marker.get()->_spillType = _spillType;

      double period = 1/_systemClockSpeed*1000;

      if (_spillType == EventWindowMarker::SpillType::offspill){
        pbtmc.get()->pbtime_ = 0;
        pbt.get()->pbtime_ = 0;
        pbt.get()->pbterr_ = 0;
        marker.get()->_eventLength = _offSpillLength*period;
      }else{
        // calculate which bin we are in from event number
        int bin = event.id().event() % _onSpillBins;
        // determine phase for this event
        // each microbunch the proton bunch gets shifted back 5 ns from the event window
        double shiftPer = -1*period/_onSpillBins;
        double phaseShift = (_initialPhaseShift + bin*shiftPer);
        if (phaseShift < -1*period)
          phaseShift += period;
        // when phase shift would be > than one clock period by next microbunch
        // the event window is one clock tick shorter
        if (phaseShift < -1*period - shiftPer)
          marker.get()->_eventLength = (_onSpillMaxLength-1)*period;
        else
          marker.get()->_eventLength = _onSpillMaxLength*period;

        pbtmc.get()->pbtime_ = _constantOffset + phaseShift;
        if (_recoFromMCTruth){
          pbt.get()->pbtime_ = pbtmc.get()->pbtime_;
          if (_recoFromMCTruthErr > 0)
            pbt.get()->pbtime_ += _randgauss.fire(0,_recoFromMCTruthErr); 
          pbt.get()->pbterr_ = _recoFromMCTruthErr;
        }else{
          pbt.get()->pbtime_ = -1*period/2.;
          pbt.get()->pbterr_ = period/sqrt(12);
        }
      }

      event.put(std::move(marker));
      event.put(std::move(pbtmc));
      event.put(std::move(pbt));
  }
}

using mu2e::EventWindowMarkerProducer;
DEFINE_ART_MODULE(EventWindowMarkerProducer);
