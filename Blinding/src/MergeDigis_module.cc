// Ed Callaghan
// Collate digis from multiple sources, and apply collision resolution
// August 2024

// stl
#include <string>
#include <vector>

// art
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

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
  class MergeDigis: public art::EDProducer{
    public:
      struct Config{
        fhicl::Sequence<art::InputTag> tracker_digi_tags{
          fhicl::Name("StrawDigiCollections"),
          fhicl::Comment("art::InputTags of source StrawDigi and StrawDigiADCWaveforms")
        };
        fhicl::Atom<bool> tracker_mc{
          fhicl::Name("MergeStrawDigiMCs"),
          fhicl::Comment("True/false to merge tracker MC truth")
        };
      };

      using Parameters = art::EDProducer::Table<Config>;
      MergeDigis(const Parameters&);

    protected:
      // tracker
      std::vector<art::InputTag> _tracker_digi_tags;
      bool _tracker_mc;
      ProditionsHandle<StrawElectronics> _tracker_conditions_handle;

    private:
      void produce(art::Event&);
  };

  // constructor
  MergeDigis::MergeDigis(const Parameters& config):
      art::EDProducer(config),
      _tracker_digi_tags(config().tracker_digi_tags()),
      _tracker_mc(config().tracker_mc()){
    // tracker
    for (const auto& tag: _tracker_digi_tags){
      this->consumes<StrawDigiCollection>(tag);
      this->consumes<StrawDigiADCWaveformCollection>(tag);
    }
    this->produces<StrawDigiCollection>();
    this->produces<StrawDigiADCWaveformCollection>();

    // calorimeter...

    // crv...

    // mc truth
    if (_tracker_mc){
      for (const auto& tag: _tracker_digi_tags){
        this->consumes<StrawDigiMCCollection>(tag);
      }
      this->produces<StrawDigiMCCollection>();
    }
  }

  void MergeDigis::produce(art::Event& event){
    // tracker: two easy steps:
    //   i) read all digis into a StrawDigiBundleCollection
    //  ii) defer collision resolution to that collection
    StrawDigiBundleCollection bundles;
    for (const auto& tag: _tracker_digi_tags){
      auto digi_handle = event.getValidHandle<StrawDigiCollection>(tag);
      auto adcs_handle = event.getValidHandle<StrawDigiADCWaveformCollection>(tag);
      if (_tracker_mc){
        auto dgmc_handle = event.getValidHandle<StrawDigiMCCollection>(tag);
        bundles.Append(*digi_handle, *adcs_handle, *dgmc_handle);
      }
      else{
        bundles.Append(*digi_handle, *adcs_handle);
      }
    }
    const auto& electronics = _tracker_conditions_handle.get(event.id());
    StrawDigiBundleCollection resolved;
    bundles.ResolveCollisions(electronics, resolved);
    auto digis = resolved.GetStrawDigiPtrs();
    auto adcs = resolved.GetStrawDigiADCWaveformPtrs();
    event.put(std::move(digis));
    event.put(std::move(adcs));
    if (_tracker_mc){
      auto dgmc = resolved.GetStrawDigiMCPtrs();
      event.put(std::move(dgmc));
    }

    // calorimeter...

    // crv...
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::MergeDigis);
