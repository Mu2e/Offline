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
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/CaloMC/inc/CaloDigiWrapperCollection.hh"

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
        fhicl::Sequence<art::InputTag> calo_digi_tags{
          fhicl::Name("CaloDigiCollections"),
          fhicl::Comment("art::InputTags of source CaloDigis")
        };
        fhicl::Atom<unsigned int> calo_adc_bits{
          fhicl::Name("CalorimeterADCBitDepth"),
          fhicl::Comment("Bit depth of calorimeter adc readings (temporary)")
        };
      };

      using Parameters = art::EDProducer::Table<Config>;
      MergeDigis(const Parameters&);

    protected:
      // tracker
      std::vector<art::InputTag> _tracker_digi_tags;
      bool _tracker_mc;
      ProditionsHandle<StrawElectronics> _tracker_conditions_handle;

      // calorimeter
      std::vector<art::InputTag> _calo_digi_tags;
      CaloDigiWrapper::sample_t _calo_max_adc;

    private:
      void produce(art::Event&);
  };

  // constructor
  MergeDigis::MergeDigis(const Parameters& config):
      art::EDProducer(config),
      _tracker_digi_tags(config().tracker_digi_tags()),
      _tracker_mc(config().tracker_mc()),
      _calo_digi_tags(config().calo_digi_tags()),
      _calo_max_adc(1 << config().calo_adc_bits()){
    // tracker
    for (const auto& tag: _tracker_digi_tags){
      this->consumes<StrawDigiCollection>(tag);
      this->consumes<StrawDigiADCWaveformCollection>(tag);
    }
    this->produces<StrawDigiCollection>();
    this->produces<StrawDigiADCWaveformCollection>();

    // calorimeter...
    for (const auto& tag: _calo_digi_tags){
      this->consumes<CaloDigiCollection>(tag);
    }
    this->produces<CaloDigiCollection>();

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

    // calorimeter: two easy steps:
    //   i) read all digis into a CaloDigiWrapperCollection
    //  ii) defer collision resolution to that collection
    CaloDigiWrapperCollection wrappers;
    for (const auto& tag: _calo_digi_tags){
      auto handle = event.getValidHandle<CaloDigiCollection>(tag);
      wrappers.Append(*handle);
    }
    CaloDigiWrapperCollection calo_resolved;
    wrappers.ResolveCollisions(_calo_max_adc, calo_resolved);
    auto calo_digis = calo_resolved.GetDigis();
    event.put(std::move(calo_digis));

    // crv...
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::MergeDigis);
