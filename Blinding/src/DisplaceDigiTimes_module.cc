// Ed Callaghan
// Deterministically shift times of all digis s.t. uniformly distributed start times are exponentially distributed
// August 2024

// stl
#include <algorithm>
#include <string>

// art
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/detail/EngineCreator.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// cetlib_except
#include "cetlib_except/exception.h"

// clhep
#include "CLHEP/Random/RandGeneral.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"

// mu2e
#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/TrackerConditions/inc/StrawElectronics.hh"
#include "Offline/TrackerMC/inc/StrawDigiBundle.hh"
#include "Offline/TrackerMC/inc/StrawDigiBundleCollection.hh"

namespace mu2e{
  class DisplaceDigiTimes: public art::EDProducer{
    public:
      struct Config{
        fhicl::Atom<art::InputTag> tracker_digi_tag{
          fhicl::Name("StrawDigiCollection"),
          fhicl::Comment("art::InputTag of digis to displace")
        };
        fhicl::DelegatedParameter target_distribution{
          fhicl::Name("path"),
          fhicl::Comment("Path to target timing distribution tabulation")
        };
      };

      using Parameters = art::EDProducer::Table<Config>;
      DisplaceDigiTimes(const Parameters&);

      void produce(art::Event&);

    protected:
      art::InputTag _tracker_digi_tag;
      ProditionsHandle<StrawElectronics> _tracker_conditions_handle;
      art::RandomNumberGenerator::base_engine_t& _engine;
      BinnedSpectrum _target_distribution;
      std::unique_ptr<CLHEP::RandGeneral> _target_sampler;
      double sample_target_time();

      TrkTypes::TDCValue straw_analog_to_digital(double, const StrawElectronics&);
      double straw_digital_to_analog(TrkTypes::TDCValue, const StrawElectronics&);

    private:
      /**/
  };

  DisplaceDigiTimes::DisplaceDigiTimes(const Parameters& config):
      art::EDProducer(config),
      _tracker_digi_tag(config().tracker_digi_tag()),
      _engine{createEngine(art::ServiceHandle<SeedService>()->getSeed())}{
    // initialize target distribution
    const auto pset = config().target_distribution.get<fhicl::ParameterSet>();
    _target_distribution = BinnedSpectrum(pset);
    _target_sampler = std::make_unique<CLHEP::RandGeneral>(_engine,
                                               _target_distribution.getPDF(),
                                               _target_distribution.getNbins());

    // tracker
    this->consumes<StrawDigiCollection>(_tracker_digi_tag);
    this->consumes<StrawDigiADCWaveformCollection>(_tracker_digi_tag);
    this->produces<StrawDigiCollection>();
    this->produces<StrawDigiADCWaveformCollection>();

    // calorimeter...
  }

  void DisplaceDigiTimes::produce(art::Event& event){
    // tracker
    auto digi_handle = event.getValidHandle<StrawDigiCollection>(_tracker_digi_tag);
    auto adcs_handle = event.getValidHandle<StrawDigiADCWaveformCollection>(_tracker_digi_tag);
    if (digi_handle->size() != adcs_handle->size()){
      std::string msg = "mismatched StrawDigiCollection and StrawDigiWaveformCollection.";
      msg += " sizes: " + std::to_string(digi_handle->size())
           + " and " + std::to_string(adcs_handle->size());
      throw cet::exception("DisplaceDigiTimes") << msg << std::endl;
    }

    // TODO currently taking first straw digi time as reference
    // this should be improved (hell, this isn't even correct!)
    // to the time of the muon decay, i.e. the time of the intersection
    // of the track with the stopping target
    const auto& electronics = _tracker_conditions_handle.get(event.id());
    auto digis = std::make_unique<StrawDigiCollection>();
    auto adcss = std::make_unique<StrawDigiADCWaveformCollection>();

    double shift = 0.0;
    if (0 < digi_handle->size()){
      auto first = digi_handle->at(0).TDC();
      auto minit = std::min_element(first.begin(), first.end());
      auto earliest_tdc = static_cast<TrkTypes::TDCValue>(*minit);
      double reference_time = this->straw_digital_to_analog(earliest_tdc, electronics);
      double target_time = this->sample_target_time();
      shift = target_time - reference_time;
    }

    // tracker
    if (0 < digi_handle->size()){
      for (const auto& old: *digi_handle){
        auto sid = old.strawId();
        auto tdc = old.TDC();
        auto tot = old.TOT();
        auto pmp = old.PMP();

        // displace tdc times
        for (size_t i = 0 ; i < TrkTypes::NENDS ; i++){
          auto analog = this->straw_digital_to_analog(tdc[i], electronics);
          auto shifted = analog + shift;
          tdc[i] = this->straw_analog_to_digital(shifted, electronics);
        }

        // new digi is identical to input digi, with displaced tdc
        digis->emplace_back(sid, tdc, tot, pmp);
      }

      // TODO indexes should be adjusted here, so that the waveforms
      // exist at the times corresponding to the tdc values
      for (const auto& old: *adcs_handle){
        adcss->push_back(old);
      }
    }
    event.put(std::move(digis));
    event.put(std::move(adcss));

    // calorimeter...
  }

  double DisplaceDigiTimes::sample_target_time(){
    const double scaled = _target_sampler->fire();
    const double rv = _target_distribution.sample(scaled);
    return rv;
  }

  TrkTypes::TDCValue DisplaceDigiTimes::straw_analog_to_digital(double time, const StrawElectronics& electronics){
    auto rv = electronics.tdcResponse(time);
    return rv;
  }

  double DisplaceDigiTimes::straw_digital_to_analog(TrkTypes::TDCValue time, const StrawElectronics& electronics){
    auto rv = static_cast<double>(time) * electronics.tdcLSB();
    return rv;
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::DisplaceDigiTimes)
