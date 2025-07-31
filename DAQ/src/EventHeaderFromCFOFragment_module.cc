// ======================================================================
//
// Add EventHeader to the event from CFO_Events
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/CFOEventFragment.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include <artdaq-core-mu2e/Data/EventHeader.hh>

#include <artdaq-core/Data/Fragment.hh>

#include <iostream>

#include <string>

#include <array>
#include <list>
#include <memory>
#include <unordered_map>
#include <vector>

#include "trace.h"
#define TRACE_NAME "Mu2eSubEventReceiver"

namespace art {
class EventHeaderFromCFOFragment;
}

// ======================================================================

class art::EventHeaderFromCFOFragment : public EDProducer {

public:
  struct Config {
    fhicl::Atom<art::InputTag> cfoTag{fhicl::Name("cfoTag"), fhicl::Comment("Input module")};
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level")};
  };

  // --- C'tor/d'tor:
  explicit EventHeaderFromCFOFragment(const art::EDProducer::Table<Config>& config);
  ~EventHeaderFromCFOFragment() override {}

  void beginRun(art::Run&) override;

  // --- Production:
  void produce(Event&) override;

private:
  art::InputTag cfoFragmentTag_;
  int diagLevel_;
};

// ======================================================================

void art::EventHeaderFromCFOFragment::beginRun(art::Run& Run) {}

art::EventHeaderFromCFOFragment::EventHeaderFromCFOFragment(
    const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    cfoFragmentTag_(config().cfoTag()), diagLevel_(config().diagLevel()) {
  produces<mu2e::EventHeader>();
}

// ----------------------------------------------------------------------

void art::EventHeaderFromCFOFragment::produce(Event& event) {

  // Collection of CaloHits for the event
  std::unique_ptr<mu2e::EventHeader> evtHdr(new mu2e::EventHeader);
  art::Handle<artdaq::Fragments> cfoFragmentHandle;

  if (!event.getByLabel(cfoFragmentTag_, cfoFragmentHandle)) {
    event.put(std::move(evtHdr));
    TLOG(TLVL_DEBUG) << "No CFO fragments found";
    return;
  }

  const auto* fragments = cfoFragmentHandle.product();
  if (fragments->size() > 0) {
    const auto& frag = fragments->at(0);
    mu2e::CFOEventFragment cfoFrag(frag);
    const CFOLib::CFO_Event cfo = cfoFrag.getData();
    // const CFO_EventRecord&   cfoRecord = cfo.GetEventRecord();
    evtHdr->mode = 0; // cfo.GetEventMode();
    evtHdr->ewt =
        static_cast<long unsigned int>(cfo.GetEventWindowTag().GetEventWindowTag().to_ullong());
    evtHdr->flags = cfo.GetEventMode().isOnSpillFlagSet();
    TLOG(TLVL_DEBUG + 20) << "mode = " << evtHdr->mode << " ewt  = " << evtHdr->ewt
                          << " flags = " << evtHdr->flags;
  } else {
    TLOG(TLVL_DEBUG) << "No CFO fragments found in the event";
  }

  event.put(std::move(evtHdr));
} // produce()

// ======================================================================

DEFINE_ART_MODULE(art::EventHeaderFromCFOFragment)

// ======================================================================
