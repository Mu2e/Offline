// ======================================================================
//
// STMPrintFragments_plugin:  Add CRV data products to the event
//
// ======================================================================

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/STMFragment.hh"
#include <artdaq-core/Data/Fragment.hh>

#include <iostream>
#include <string>
#include <memory>

namespace art
{
  class STMPrintFragments;
}

using art::STMPrintFragments;

// ======================================================================

class art::STMPrintFragments : public EDAnalyzer
{
  public:
  struct Config
   {
  //   //    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic Level")};
  //   //    fhicl::Atom<art::InputTag> CRVDataDecodersTag{fhicl::Name("crvTag"),
  //   //                                               fhicl::Comment("crv Fragments Tag")};
  };

  // --- C'tor/d'tor:
  explicit STMPrintFragments(const art::EDAnalyzer::Table<Config>& config);

  // --- Production:
  virtual void analyze(const Event&);

  private:
  //  int decompressCrvDigi(uint8_t adc);
  //  int16_t decompressCrvDigi(int16_t adc);

  //  int                                      _diagLevel;
  //  art::InputTag                            _CRVDataDecodersTag;
  //  mu2e::ProditionsHandle<mu2e::CRVOrdinal> _channelMap_h;

}; // STMPrintFragments

// ======================================================================

STMPrintFragments::STMPrintFragments(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config}
{
  //  produces<mu2e::CrvDigiCollection>();
}

// ----------------------------------------------------------------------

void STMPrintFragments::analyze(const Event& event)
{
  art::EventNumber_t eventNumber = event.event();

  auto STMFragments = event.getValidHandle<artdaq::Fragments>("daq:STM");

  std::cout << std::dec << "Analyzer: Run " << event.run() << ", subrun " << event.subRun()
            << ", event " << eventNumber << " has " << std::endl;
  std::cout << STMFragments->size() << " STM fragments." << std::endl;

  // std::vector<art::Handle<std::vector<artdaq::Fragment>>> fragmentHandles;
  // fragmentHandles = event.getMany<std::vector<artdaq::Fragment>>();
  // std::cout << "AE: fragmentHandles.size() = " << fragmentHandles.size() << std::endl;
  // for (auto const& hndl : fragmentHandles) {
  //   std::cout << "AE: hdnle = " << hndl << std::endl;
  //   for (auto const& fragment : *hndl) {
  //       int fragID = fragment.fragmentID();
  //       std::cout << "AE: fragID = " << fragID << std::endl;
  //   }
  // }
  int frag_counter = 0;
  for (auto& frag : *STMFragments) {
    ++frag_counter;

    //    auto stm_frag = static_cast<mu2e::STMFragment>(frag);
    const auto dataBegin = frag.dataBegin();
    const auto dataEnd = frag.dataEnd();
    const auto stmDataBegin = reinterpret_cast<int16_t const*>(dataBegin);
    const auto stmDataEnd = reinterpret_cast<int16_t const*>(dataEnd);
    auto frag_id = frag.fragmentID();
    std::cout << "frag_id = " << frag_id << std::endl;
    for (auto i = stmDataBegin; i != stmDataEnd; ++i) {
      std::cout << "Frag #" << frag_counter << ": *(stmDataBegin+" << i - stmDataBegin << ") = " << *i << std::endl;
    }

    // //    std::cout << "Trigger Header Address: " << stm_frag.GetTHdr() << std::endl;
    // std::cout << "Frag #" << frag_counter << ": EvNum: " << *(stm_frag.EvNum()) << std::endl;
    // std::cout << "Frag #" << frag_counter << ": DataType: " << *(stm_frag.DataType()) << std::endl;
    // std::cout << "Frag #" << frag_counter << ": EvLen: " << *(stm_frag.EvLen()) << std::endl;

    // unsigned int max_samples = (*(stm_frag.EvLen()))/100.;
    // std::cout << "Frag #" << frag_counter << ": First " << max_samples << " int16s of data: ";
    // for (size_t i = 0; i < max_samples; ++i) {
    //   std::cout << *(stm_frag.DataBegin()+i) << " ";
    // }
    // std::cout << std::endl;
    //    std::cout << "Trigger Header Channel: " << *(stm_frag.GetTHdr()) << std::endl;
    //    std::cout << "Trigger Header EvNum: " << *(stm_frag.GetTHdr()+8) << std::endl;
  }

} // produce()

// ======================================================================

DEFINE_ART_MODULE(STMPrintFragments)
