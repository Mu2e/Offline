// ======================================================================
//
// STMDigisFromFragments: create all types of STMDigis from STMFragments
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/RecoDataProducts/inc/STMMWDDigi.hh"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/STMFragment.hh"
#include <artdaq-core/Data/Fragment.hh>

#include <string>
#include <memory>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <list>
#include <numeric>
#include <random>
#include <vector>

// ROOT includes
#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TString.h>
#include <TLine.h>
#include <TGraph.h>

namespace art
{
  class STMDigisFromFragments;
}

using art::STMDigisFromFragments;

// ======================================================================

class art::STMDigisFromFragments : public EDProducer
{
  public:
  struct Config
  {
    fhicl::Atom<art::InputTag> stmTag {fhicl::Name("stmTag"), fhicl::Comment("Input module")};
    // TODO: add fhicl parameters to that we can choose which types of fragments we read out
    // e.g.     fhicl::Atom<bool> processRaw {fhicl::Name("processRaw"), fhicl::Comment("Process Raw STMFragments")};
  };

  // --- C'tor/d'tor:
  explicit STMDigisFromFragments(const art::EDProducer::Table<Config>& config);

  // --- Production:
  virtual void produce(Event&);

  private:

  art::InputTag _stmFragmentsTag;

}; // STMDigisFromFragments

// ======================================================================


STMDigisFromFragments::STMDigisFromFragments(const art::EDProducer::Table<Config>& config) :
  art::EDProducer{config}
  ,_stmFragmentsTag(config().stmTag())

{
  // Set the size of the vector
  produces<mu2e::STMWaveformDigiCollection>("raw");
  produces<mu2e::STMWaveformDigiCollection>("zs");
  produces<mu2e::STMWaveformDigiCollection>("mwd"); // TODO: we should create an STMMWDDigi collection instead of STMWaveformDigis for this
}

// ----------------------------------------------------------------------


void STMDigisFromFragments::produce(Event& event)
{
  std::unique_ptr<mu2e::STMWaveformDigiCollection> raw_waveform_digis(new mu2e::STMWaveformDigiCollection);
  std::unique_ptr<mu2e::STMWaveformDigiCollection> zs_waveform_digis(new mu2e::STMWaveformDigiCollection);
  std::unique_ptr<mu2e::STMWaveformDigiCollection> mwd_waveform_digis(new mu2e::STMWaveformDigiCollection);

  art::Handle<artdaq::Fragments> STMFragmentsH;
  event.getByLabel(_stmFragmentsTag, STMFragmentsH);
  const auto STMFragments = STMFragmentsH.product();

  for (const auto& frag : *STMFragments) {
    const auto& stm_frag = static_cast<mu2e::STMFragment>(frag);

    mu2e::STMWaveformDigi stm_waveform;

    if (stm_frag.isRaw()) {
      stm_waveform.set_data(stm_frag.payloadWords(), stm_frag.payloadBegin());
      raw_waveform_digis->emplace_back(stm_waveform);
    }
    else if (stm_frag.isZS()) {
      stm_waveform.set_data(stm_frag.payloadWords(), stm_frag.payloadBegin());
      zs_waveform_digis->emplace_back(stm_waveform);
    }
    else if (stm_frag.isMWD()) {
      stm_waveform.set_data(stm_frag.payloadWords(), stm_frag.payloadBegin());
      mwd_waveform_digis->emplace_back(stm_waveform);
    }
  }

  event.put(std::move(raw_waveform_digis), "raw");
  event.put(std::move(zs_waveform_digis), "zs");
  event.put(std::move(mwd_waveform_digis), "mwd");

} // produce()

// ======================================================================

DEFINE_ART_MODULE(STMDigisFromFragments)
