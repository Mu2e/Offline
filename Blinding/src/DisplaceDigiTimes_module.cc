// Ed Callaghan
// Coherently shift times of all digis s.t. event timing follows a tabulated distribution
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
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Offline/RecoDataProducts/inc/KalIntersection.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
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
        fhicl::Atom<art::InputTag> kalseed_tag{
          fhicl::Name("KalSeedCollection"),
          fhicl::Comment("art::InputTag of KalSeedCollection to prescale")
        };
        fhicl::DelegatedParameter target_distribution{
          fhicl::Name("target_distribution"),
          fhicl::Comment("BinnedSpectrum configuration")
        };
        fhicl::Sequence<std::string> surface_ids{
          fhicl::Name("SurfaceIds"),
          fhicl::Comment("Prioritized sequence of mu2e::SurfaceIds at which tracks may be sampled for reference timing")
        };
      };

      using Parameters = art::EDProducer::Table<Config>;
      DisplaceDigiTimes(const Parameters&);

      void produce(art::Event&);

    protected:
      art::InputTag _tracker_digi_tag;
      art::InputTag _kalseed_tag;
      SurfaceIdCollection _surface_ids;
      ProditionsHandle<StrawElectronics> _tracker_conditions_handle;
      art::RandomNumberGenerator::base_engine_t& _engine;
      BinnedSpectrum _target_distribution;
      std::unique_ptr<CLHEP::RandGeneral> _target_sampler;
      double sample_target_time();

    private:
      /**/
  };

  DisplaceDigiTimes::DisplaceDigiTimes(const Parameters& config):
      art::EDProducer(config),
      _tracker_digi_tag(config().tracker_digi_tag()),
      _kalseed_tag(config().kalseed_tag()),
      _engine{createEngine(art::ServiceHandle<SeedService>()->getSeed())}{
    // initialize target distribution
    const auto pset = config().target_distribution.get<fhicl::ParameterSet>();
    _target_distribution = BinnedSpectrum(pset);
    _target_sampler = std::make_unique<CLHEP::RandGeneral>(_engine,
                                               _target_distribution.getPDF(),
                                               _target_distribution.getNbins());

    // surfaces at which timing may be referenced
    for (const auto& name: config().surface_ids()){
      _surface_ids.emplace_back(name);
    }

    // reco
    this->consumes<KalSeed>(_kalseed_tag);

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

    // first, check that there is only a single KalSeed associated with
    // these digis
    auto kalseeds_handle = event.getValidHandle<KalSeedCollection>(_kalseed_tag);
    if (kalseeds_handle->size() != 1){
      std::string msg = "ambiguous KalSeedCollection (length " + std::to_string(kalseeds_handle->size()) + "), cannot assign reference time";
      throw cet::exception("DisplaceDigiTimes") << msg << std::endl;
    }
    auto kalseed = (*kalseeds_handle)[0];

    // search for downward-going KalIntersection with a matching surface
    KalIntersection intersection;
    auto sit = _surface_ids.begin();
    bool adequate = false;
    while ((!adequate) && (sit != _surface_ids.end())){
      const auto& iits = kalseed.intersections(*sit);
      // there are at most 2 valid intersections (upgoing / downgoing)
      // so selecting on which one is downgoing is sufficient
      for (const auto& iit: iits){
        auto direction = iit->momentum3().Unit().Z();
        if (0 < direction){
          intersection = *iit;
          adequate = true;
        }
      }
      sit++;
    }

    // if no valid intersection found, then we cannot use this event
    if (!adequate){
      std::string msg = "no valid KalIntersection found";
      throw cet::exception("DisplaceDigiTimes") << msg << std::endl;
    }
    // otherwise, take the time of the intersection as the reference time
    double reference_time = intersection.time();
    double target_time = this->sample_target_time();
    double shift = target_time - reference_time;

    // apply the shift to tracker digis
    const auto& electronics = _tracker_conditions_handle.get(event.id());
    auto digis = std::make_unique<StrawDigiCollection>();
    auto adcss = std::make_unique<StrawDigiADCWaveformCollection>();
    if (0 < digi_handle->size()){
      for (const auto& old: *digi_handle){
        auto sid = old.strawId();
        auto tdc = old.TDC();
        auto tot = old.TOT();
        auto pmp = old.PMP();

        // displace tdc times
        for (size_t i = 0 ; i < TrkTypes::NENDS ; i++){
          auto analog = electronics.timeDigitalToAnalog(tdc[i], sid);
          auto shifted = analog + shift;
          tdc[i] = electronics.timeAnalogToDigital(shifted, sid);
        }

        // new digi is identical to input digi, with displaced tdc
        digis->emplace_back(sid, tdc, tot, pmp);
      }

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
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::DisplaceDigiTimes)
