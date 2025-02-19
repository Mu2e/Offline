// Ed Callaghan
// Generate uncorrelated tracker digis assuming a poisson rate of MIP pulses
// February 2025

// stl
#include <memory>

// art
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"

// canvas
#include "canvas/Utilities/InputTag.h"

// clhep
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

// cetlib_except
#include "cetlib_except/exception.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/TrackerConditions/inc/StrawElectronics.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerMC/inc/AnalogWireSignal.hh"
#include "Offline/TrackerMC/inc/AnalogSignalShapeTool.hh"

namespace mu2e{
  class PoissonTrackerNoise: public art::EDProducer{
    public:
      struct Config{
        fhicl::Atom<double> rate{
          fhicl::Name("rate"),
          fhicl::Comment("Per-straw rate (in GHz) of pulses")
        };
        fhicl::DelegatedParameter signal{
          fhicl::Name("signal"),
          fhicl::Comment("Tool configuration to produce primary signal shapes");
        };
        fhicl::Atom<art::InputTag> ewm_tag{
          fhicl::Name("EventWindowMarker"),
          fhicl::Comment("art::InputTag of EventWindowMarker for event")
        };
        fhicl::Atom<double> crtime_spacing{
          fhicl::Name("threshold_crossing_coarse_sampling"),
          fhicl::Comment("Resolution of coarse linear scan for threshold crossing")
        };
        fhicl::Atom<double> crtime_tolerance{
          fhicl::Name("threshold_crossing_tolerance"),
          fhicl::Comment("Tolerance on exact threshold-crossing time")
        };
      };

      using Parameters = art::EDProducer::Table<Config>;
      PoissonTrackerNoise(const Parameters&);
     ~PoissonTrackerNoise() = default;

      void produce(art::Event&) override;

    protected:
      double _rate;
      art::InputTag _ewm_tag;
      double _crtime_spacing;
      double _crtime_tolerance;
      ProditionsHandle<StrawElectronics> _electronics_h;
      ProditionsHandle<StrawResponse> _response_h;
      ProditionsHandle<Tracker> _tracker_h;
      art::RandomNumberGenerator::base_engine_t& _engine;
      std::unique_ptr<CLHEP::RandPoisson> _poisson;       // for normalization
      std::unique_ptr<CLHEP::RandFlat> _uniform;          // for timing
      std::unique_ptr<AnalogSignalShapeTool> _shape;      // signal shape

    private:
      /**/
  };

  PoissonTrackerNoise::PoissonTrackerNoise(const Parameters& config):
      art::EDProducer(config),
      _rate(config().rate()),
      _ewm_tag(config().ewm_tag()),
      _crtime_spacing(config().crtime_spacing()),
      _crtime_tolerance(config().crtime_tolerance()),
      _engine{createEngine(art::ServiceHandle<SeedService>()->getSeed())}{

    // rng initializations
    _uniform = std::make_unique<CLHEP::RandFlat>(_engine);
    _poisson = std::make_unique<CLHEP::RandPoisson>(_engine);

    // signal-shape sampler
    auto signal_config = config().signal.get<fhicl::ParameterSet>();
    _shape = art::make_tool<AnalogSignalShapeTool>(signal_config);

    // framework hooks
    this->produces<StrawDigiCollection>();
    this->produces<StrawDigiADCWaveformCollection>();
    this->produces<StrawDigiMCCollection>();
  }

  void PoissonTrackerNoise::produce(art::Event& event){
    // fetch straw conditions
    const Tracker& tracker = _tracker_h.get(event.id());
    const StrawResponse& response = _response_h.get(event.id());
    const StrawElectronics& electronics = _electronics_h.get(event.id());
    const auto internal_delay = electronics.electronicsTimeDelay();

    // define analog event window, beginning at time == 0
    auto ewm_h = event.getValidHandle<EventWindowMarker>(_ewm_tag);
    double window = ewm_h->eventLength();
    TrkTypes::TDCValue maxTDC = electronics.tdcResponse(window - internal_delay);

    // used downstream
    double threshold;

    // containers for payload
    auto digis = std::make_unique<StrawDigiCollection>();
    auto adcss = std::make_unique<StrawDigiADCWaveformCollection>();
    auto dgmcs = std::make_unique<StrawDigiMCCollection>();

    // loop over all straws
    for (uint16_t iplane = 0 ; iplane < StrawId::_nplanes ; iplane++){
      for (uint16_t ipanel = 0 ; ipanel < StrawId::_npanels ; ipanel++){
        for (uint16_t istraw = 0 ; istraw < StrawId::_nstraws ; istraw++){
          // TODO tracker status here?
          auto sid = StrawId(iplane, ipanel, istraw);
          const Straw& straw = tracker.getStraw(sid);

          // sample number of pulses that this straw will see
          double mean = _rate * window;
          size_t count = static_cast<size_t>(_poisson->fire(mean));

          // loop over all pulses
          for (size_t i = 0 ; i < count ; i++){
            // first, assign time that signal is induced on wire
            double nominal_time = _uniform->fire(0.0, window);

            // next, assign position along the wire where signal is induced...
            double length = 2.0 * straw.halfLength();
            double position = _uniform->fire(0.0, length);

            // ...which translates into two times at either straw end
            // TODO units of prop speed mm /ns -> cm / ns ?
            // TODO nominal energy scale here...?
            double edep = 0.5 * 5.9; // half of 55Fe 5.9 keV line
            double propagation_speed = 2.0 * response.halfPropV(sid, edep);
            double lht = nominal_time +           position /propagation_speed;
            double rht = nominal_time + (length - position)/propagation_speed;

            // sample signal shape
            auto shape = _shape->Sample();
            auto signal = AnalogWireSignal(shape);

            // fill two-sided waveforms with delayed analog signal
            // TODO: transfer function to affect asymmetric transmission
            // TODO: high-frequency noise
            auto lhs = signal;
            threshold = electronics.threshold(sid, StrawEnd::cal);
            bool lct = lhs.TranslateToThresholdCrossingTime(threshold, lht,
                                                            0.0, window,
                                                            _crtime_spacing,
                                                            _crtime_tolerance);
            auto rhs = signal;
            threshold = electronics.threshold(sid, StrawEnd::hv);
            bool rct = rhs.TranslateToThresholdCrossingTime(threshold, rht,
                                                            0.0, window,
                                                            _crtime_spacing,
                                                            _crtime_tolerance);

            // if both sides are above threshold, try to digitize
            bool two_sided = lct && rct;
            if (two_sided){
              // forward-declare downstream digital vessels
              TrkTypes::TDCTimes crTimesPhysical;     // threshold-crossing times
              TrkTypes::TDCValues crTimesDigital;     // ...after digitizing
            //TrkTypes::TOTTimes atTimesPhysical;     // times over threshold
              TrkTypes::TOTValues atTimesDigital;     // ...after digitizing
              TrkTypes::ADCTimes wfTimes;             // times of wf samples
            //TrkTypes::ADCVoltages samplesPhysical;  // analog sample values
              TrkTypes::ADCWaveform samplesDigital;   // digital sample values
              TrkTypes::ADCValue pmp;                 // digital amplitude

              // absolute times must be uncalibrated and then fed through tdcs
              bool isOnSpill = (ewm_h->spillType() == EventWindowMarker::onspill);
              // by definition, from the above
              crTimesPhysical[StrawEnd::cal] = lht;
              crTimesPhysical[StrawEnd::hv]  = rht;
              electronics.uncalibrateTimes(crTimesPhysical, sid);
              bool in_window = electronics.digitizeTimes(crTimesPhysical,
                                                         crTimesDigital,
                                                         isOnSpill,
                                                         maxTDC);
              // if both tdcs are within the digitization window,
              // then sum waveforms and create a digi
              if (in_window){
                // calculate double-sided tots
                threshold = electronics.threshold(sid, StrawEnd::cal);
                lhs.DigitalTimeOverThreshold(electronics, threshold, lht,
                                             atTimesDigital[StrawEnd::cal]);
                threshold = electronics.threshold(sid, StrawEnd::hv);
                rhs.DigitalTimeOverThreshold(electronics, threshold, rht,
                                             atTimesDigital[StrawEnd::hv]);
                // sum analog signals
                auto summed = lhs + rhs;
                // sample digital waveform and compute amplitude
                summed.Digitize(electronics, sid, crTimesPhysical[0],
                                wfTimes, samplesDigital, pmp);
                // emplace new digi, waveform, and mc truth
                digis->emplace_back(sid, crTimesDigital, atTimesDigital, pmp);
                adcss->emplace_back(samplesDigital);
                dgmcs->emplace_back(sid, true);
              }
            }
          }
        }
      }
    }

    // store in event
    event.put(std::move(digis));
    event.put(std::move(adcss));
    event.put(std::move(dgmcs));
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::PoissonTrackerNoise)
