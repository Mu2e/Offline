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
    fhicl::Atom<art::InputTag> stmTag {fhicl::Name("stmTag"), fhicl::Comment("Input module")};
  };

  // --- C'tor/d'tor:
  explicit STMWaveformDigisFromFragments(const art::EDProducer::Table<Config>& config);

  // --- Production:
  virtual void produce(Event&);

  private:

  int eventCounter = 0;
  
  art::InputTag _stmFragmentsTag;
  
}; // STMWaveformDigisFromFragments

// ======================================================================

STMWaveformDigisFromFragments::STMWaveformDigisFromFragments(const art::EDProducer::Table<Config>& config) :
  art::EDProducer{config},
  _stmFragmentsTag(config().stmTag())
{
  produces<mu2e::STMWaveformDigiCollection>();
}

// ----------------------------------------------------------------------


void STMWaveformDigisFromFragments::produce(Event& event)
{
  eventCounter++;
  std::unique_ptr<mu2e::STMWaveformDigiCollection> stm_waveform_digis(new mu2e::STMWaveformDigiCollection);
  
  //auto STMFragments = event.getValidHandle<artdaq::Fragments>("daq:STM");
  art::Handle<artdaq::Fragments> STMFragmentsH;
  if(!event.getByLabel(_stmFragmentsTag, STMFragmentsH)) {
    event.put(std::move(stm_waveform_digis));
    return;
  }
  const auto STMFragments = STMFragmentsH.product();

  int16_t detID = 0;
  for (const auto& frag : *STMFragments) {
    const auto& stm_frag = static_cast<mu2e::STMFragment>(frag);
    //int16_t event_number = *(stm_frag.EvNum());
    int16_t event_len = *(stm_frag.EvLen());    
    detID = *(stm_frag.detID()) & 0xFF;
    //std::cout << "DetID: " << *(stm_frag.detID()) << "  |   " << detID << "\n";
    //std::cout << "In event: " << event_number << " | Length of event: " << event_len << "\n";
    std::vector<int16_t> adcs;
    adcs.reserve(event_len);
    for (unsigned long int i_adc_sample = 0; i_adc_sample < event_len; ++i_adc_sample) {
      adcs.emplace_back(*(stm_frag.DataBegin()+i_adc_sample));
    }

    // Create the STMWaveformDigi and put it in the event
    uint32_t trig_time_offset = 0;
    mu2e::STMWaveformDigi stm_waveform(detID, trig_time_offset, adcs);
    stm_waveform_digis->push_back(stm_waveform);
  }


  event.put(std::move(stm_waveform_digis));

} // produce()

// ======================================================================

DEFINE_ART_MODULE(STMWaveformDigisFromFragments)
