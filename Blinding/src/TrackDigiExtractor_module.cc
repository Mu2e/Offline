// Ed Callaghan
// Extract digis which participate in a prefit track
// September 2024

// stl
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

// mu2e
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"

namespace mu2e{
  class TrackDigiExtractor: public art::EDProducer{
    public:
      struct Config{
        // track
        fhicl::Atom<art::InputTag> kalseed_tag{
          fhicl::Name("KalSeedCollection"),
          fhicl::Comment("art::InputTag of KalSeeds to extract digis from")
        };

        // tracker
        fhicl::Atom<art::InputTag> strawhit_tag{
          fhicl::Name("StrawHitCollection"),
          fhicl::Comment("art::InputTag of StrawHits associated with KalSeed")
        };
        fhicl::Atom<art::InputTag> strawdigi_tag{
          fhicl::Name("StrawDigiCollection"),
          fhicl::Comment("art::InputTag of StrawDigis associated with StrawHits")
        };

        // calorimeter...
      };

      using Parameters = art::EDProducer::Table<Config>;
      TrackDigiExtractor(const Parameters&);

    protected:
      art::InputTag _kalseed_tag;
      art::InputTag _strawhit_tag;
      art::InputTag _strawdigi_tag;

    private:
      void produce(art::Event&);
  };

  TrackDigiExtractor::TrackDigiExtractor(const Parameters& config):
      art::EDProducer(config),
      _kalseed_tag(config().kalseed_tag()),
      _strawhit_tag(config().strawhit_tag()),
      _strawdigi_tag(config().strawdigi_tag()){
    this->consumes<KalSeedCollection>(_kalseed_tag);
    this->consumes<StrawHitCollection>(_strawhit_tag);
    this->consumes<StrawDigiCollection>(_strawdigi_tag);
    this->consumes<StrawDigiADCWaveformCollection>(_strawdigi_tag);
    this->produces<StrawDigiCollection>();
    this->produces<StrawDigiADCWaveformCollection>();
  }

  void TrackDigiExtractor::produce(art::Event& event){
    // data-product lookups
    auto kalseeds = event.getValidHandle<KalSeedCollection>(_kalseed_tag);
    auto strawhits = event.getValidHandle<StrawHitCollection>(_strawhit_tag);
    auto strawdigis = event.getValidHandle<StrawDigiCollection>(_strawdigi_tag);
    auto strawadcss = event.getValidHandle<StrawDigiADCWaveformCollection>(_strawdigi_tag);

    // sanity checks
    if (strawhits->size() != strawdigis->size()){
      std::string msg = "mismatched StrawHit and StrawDigi collection lengths";
      throw cet::exception("TrackDigiExtractor") << msg << std::endl;
    }
    if (strawdigis->size() != strawadcss->size()){
      std::string msg = "mismatched StrawDigi and StrawDigiADCWaveform collection lengths";
      throw cet::exception("TrackDigiExtractor") << msg << std::endl;
    }

    // deep-copy digis underlying tracks
    auto ext_strawdigis = std::make_unique<StrawDigiCollection>();
    auto ext_strawadcss = std::make_unique<StrawDigiADCWaveformCollection>();
    for (const auto& kalseed: *kalseeds){
      // tracker
      auto seeds = kalseed.hits();
      for (const auto& seed: seeds){
        auto idx = seed.index();
        auto hit = (*strawhits)[idx];
        auto digi = (*strawdigis)[idx];
        auto adcs = (*strawadcss)[idx];

        // sanity check
        if (hit.strawId() != digi.strawId()){
          std::string msg = "mismatched StrawHit and StrawDigi StrawIDs";
          throw cet::exception("TrackDigiExtractor") << msg << std::endl;
        }

        ext_strawdigis->push_back(digi);
        ext_strawadcss->push_back(adcs);
      }

      // calorimeter...
    }

    // place into event
    event.put(std::move(ext_strawdigis));
    event.put(std::move(ext_strawadcss));
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::TrackDigiExtractor)
