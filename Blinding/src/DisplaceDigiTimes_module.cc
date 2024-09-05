// Ed Callaghan
// Deterministically shift times of all digis s.t. uniformly distributed start times are exponentially distributed
// August 2024

// stl
#include <algorithm>
#include <string>

// art
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

// cetlib_except
#include "cetlib_except/exception.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"

// mu2e
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
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
        fhicl::Atom<double> offspill_min{
          fhicl::Name("OffSpillTimeMin"),
          fhicl::Comment("Smallest accepted OffSpill time")
        };
        fhicl::Atom<double> offspill_max{
          fhicl::Name("OffSpillTimeMax"),
          fhicl::Comment("Largest accepted OffSpill time")
        };
        fhicl::Atom<double> onspill_min{
          fhicl::Name("OnSpillTimeMin"),
          fhicl::Comment("Smallest mapped OnSpill time")
        };
        fhicl::Atom<double> onspill_max{
          fhicl::Name("OnSpillTimeMax"),
          fhicl::Comment("Largest mapped OnSpill time")
        };
        fhicl::Atom<double> onspill_lifetime{
          fhicl::Name("OnSpillLifetime"),
          fhicl::Comment("Decay time of mapped OnSpill times")
        };
      };

      using Parameters = art::EDProducer::Table<Config>;
      DisplaceDigiTimes(const Parameters&);

      void produce(art::Event&);

    protected:
      art::InputTag _tracker_digi_tag;
      ProditionsHandle<StrawElectronics> _tracker_conditions_handle;
      double _offspill_min;
      double _offspill_max;
      double _onspill_min;
      double _onspill_max;
      double _onspill_lifetime;
      double displaced_time(double);

      TrkTypes::TDCValue straw_analog_to_digital(double, const StrawElectronics&);
      double straw_digital_to_analog(TrkTypes::TDCValue, const StrawElectronics&);

    private:
      /**/
  };

  DisplaceDigiTimes::DisplaceDigiTimes(const Parameters& config):
      art::EDProducer(config),
      _tracker_digi_tag(config().tracker_digi_tag()),
      _offspill_min(config().offspill_min()),
      _offspill_max(config().offspill_max()),
      _onspill_min(config().onspill_min()),
      _onspill_max(config().onspill_max()),
      _onspill_lifetime(config().onspill_lifetime()){
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
      double target_time = this->displaced_time(reference_time);
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

  double DisplaceDigiTimes::displaced_time(double time){
    double shifted = time - _offspill_min;
    double window_length = _offspill_max - _offspill_min;
    double rv = - _onspill_lifetime * log(shifted / window_length);
    rv += _onspill_min;
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
