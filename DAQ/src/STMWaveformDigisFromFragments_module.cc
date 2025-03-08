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

    // Get 40MHz DTC clock time
    uint16_t DTC0 =  static_cast<uint16_t>((*(stm_frag.GetTHdr() + 28) >> 8) & 0xF) & 0xFFFF;
    uint16_t DTC1 =  static_cast<uint16_t>(*(stm_frag.GetTHdr() + 29)) & 0xFFFF;
    uint16_t DTC2 =  static_cast<uint16_t>(*(stm_frag.GetTHdr() + 30)) & 0xFFFF;
    uint16_t DTC3 =  static_cast<uint16_t>(*(stm_frag.GetTHdr() + 31)) & 0xFFFF;
    uint64_t DTC = static_cast<uint64_t>(DTC3) << 40 | static_cast<uint64_t>(DTC2) << 24 | static_cast<uint64_t>(DTC1) << 8 | static_cast<uint64_t>(DTC0);
    //std::cout << "DTCs: " << std::hex << DTC0 << " " << DTC1 << " " << DTC2 << " " << DTC3 << " --> " << DTC << std::dec << std::endl;
    //std::cout << "DTC time = " << DTC << " = " << (DTC/40e6) << "\n";

    // Get 75MHz ADC clock time
    uint16_t ADC0 =  static_cast<uint16_t>(*(stm_frag.GetTHdr() + 4)) & 0xFFFF;
    uint16_t ADC1 =  static_cast<uint16_t>(*(stm_frag.GetTHdr() + 5)) & 0xFFFF;
    uint16_t ADC2 =  static_cast<uint16_t>(*(stm_frag.GetTHdr() + 6)) & 0xFFFF;
    uint16_t ADC3 =  static_cast<uint16_t>(*(stm_frag.GetTHdr() + 7)) & 0xFFFF;
    uint64_t ADC = static_cast<uint64_t>(ADC3) << 48 | static_cast<uint64_t>(ADC2) << 32 | static_cast<uint64_t>(ADC1) << 16 | static_cast<uint64_t>(ADC0);
    //std::cout << "ADCs: " << std::hex << ADC0 << " " << ADC1 << " " << ADC2 << " " << ADC3 << " --> " << ADC << std::dec << std::endl;
    //std::cout << "ADC time = " << ADC << " = " << (ADC/75e6) << "\n";
    
    std::vector<int16_t> adcs;
    adcs.reserve(event_len);
    for (unsigned long int i_adc_sample = 0; i_adc_sample < event_len; ++i_adc_sample) {
      adcs.emplace_back(*(stm_frag.DataBegin()+i_adc_sample));
    }

    // Create the STMWaveformDigi and put it in the event
    uint32_t trig_time_offset = 0;
    mu2e::STMWaveformDigi stm_waveform(detID, DTC, ADC, trig_time_offset, adcs);
    stm_waveform_digis->push_back(stm_waveform);
  }


  event.put(std::move(stm_waveform_digis));

} // produce()

// ======================================================================

DEFINE_ART_MODULE(STMWaveformDigisFromFragments)
