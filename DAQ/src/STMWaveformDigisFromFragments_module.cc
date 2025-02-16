// ======================================================================
//
// STMWaveformDigisFromFragments: create STMWaveformDigis from STMFragments
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/STMFragment.hh"
#include <artdaq-core/Data/Fragment.hh>

#include <iostream>
#include <string>
#include <memory>

namespace art
{
  class STMWaveformDigisFromFragments;
}

using art::STMWaveformDigisFromFragments;

// ======================================================================

class art::STMWaveformDigisFromFragments : public EDProducer
{
  public:
  struct Config
  {
    // fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic Level")};
    // fhicl::Atom<art::InputTag> CRVDataDecodersTag{fhicl::Name("crvTag"),
    //                                            fhicl::Comment("crv Fragments Tag")};
  };

  // --- C'tor/d'tor:
  explicit STMWaveformDigisFromFragments(const art::EDProducer::Table<Config>& config);

  // --- Production:
  virtual void produce(Event&);

  private:

}; // STMWaveformDigisFromFragments

// ======================================================================

STMWaveformDigisFromFragments::STMWaveformDigisFromFragments(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}
{
  produces<mu2e::STMWaveformDigiCollection>();
}

// ----------------------------------------------------------------------


void STMWaveformDigisFromFragments::produce(Event& event)
{

  auto STMFragments = event.getValidHandle<artdaq::Fragments>("daq:STM");
  std::unique_ptr<mu2e::STMWaveformDigiCollection> stm_waveform_digis(new mu2e::STMWaveformDigiCollection);

  for (const auto& frag : *STMFragments) {
    const auto& stm_frag = static_cast<mu2e::STMFragment>(frag);
    int16_t event_number = *(stm_frag.EvNum());
    int16_t event_len = *(stm_frag.EvLen());
    std::cout << "In event: " << event_number << " | Length of event: " << event_len << "\n";
    std::vector<int16_t> adcs;
    adcs.reserve(event_len);
    for (unsigned long int i_adc_sample = 0; i_adc_sample < event_len; ++i_adc_sample) {
      adcs.emplace_back(*(stm_frag.DataBegin()+i_adc_sample));
    }

    // Create the STMWaveformDigi and put it in the event
    uint32_t trig_time_offset = 0;
    mu2e::STMWaveformDigi stm_waveform(trig_time_offset, adcs);
    stm_waveform_digis->push_back(stm_waveform);
  }


  event.put(std::move(stm_waveform_digis));

} // produce()

// ======================================================================

DEFINE_ART_MODULE(STMWaveformDigisFromFragments)
