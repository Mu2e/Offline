// ======================================================================
//
// Add EventHeader to the event from CFO_Events
//
// ======================================================================

// Framework
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"

// artdaq
#include <artdaq-core-mu2e/Data/EventHeader.hh>
#include "artdaq-core-mu2e/Overlays/CFOEventFragment.hh"
#include "artdaq-core-mu2e/Overlays/CFO_Packets/CFO_Event.h"
#include "artdaq-core-mu2e/Overlays/CFO_Packets/CFO_EventRecord.h"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include <artdaq-core/Data/Fragment.hh>

// Offline
#include "Offline/DataProducts/inc/EventWindowMarker.hh"

// C++
#include <iostream>
#include <string>
#include <array>
#include <list>
#include <memory>
#include <unordered_map>
#include <vector>

// TRACE
#include "trace.h"
#define TRACE_NAME "EventHeaderFromCFOFragment"

namespace art {
class EventHeaderFromCFOFragment;
}

// ======================================================================

class art::EventHeaderFromCFOFragment : public EDProducer {

public:
  struct Config {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<art::InputTag> cfoTag   {Name("cfoTag"),    Comment("Input CFO fragment tag")};
    fhicl::Atom<bool>          ewm      {Name("createEWM"), Comment("Produce an event window marker object")};
    fhicl::Atom<int>           diagLevel{Name("diagLevel"), Comment("diagnostic level")};
  };

  // --- C'tor/d'tor:
  explicit EventHeaderFromCFOFragment(const art::EDProducer::Table<Config>& config);
  ~EventHeaderFromCFOFragment() override {}

  void beginRun(art::Run&) override;

  // --- Production:
  void produce(Event&) override;

private:
  art::InputTag cfoFragmentTag_;
  bool          ewm_;
  int           diagLevel_;
};

// ======================================================================

void art::EventHeaderFromCFOFragment::beginRun(art::Run& Run) {}

art::EventHeaderFromCFOFragment::EventHeaderFromCFOFragment(
    const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    cfoFragmentTag_(config().cfoTag()),
    ewm_(config().ewm()),
    diagLevel_(config().diagLevel()) {
  produces<mu2e::EventHeader>();
  if(ewm_) {
    TLOG(TLVL_DEBUG + 2) << "Producing EventWindowMarker";
    produces<mu2e::EventWindowMarker>();
  }
}

// ----------------------------------------------------------------------

void art::EventHeaderFromCFOFragment::produce(Event& event) {

  // Collection of CaloHits for the event
  std::unique_ptr<mu2e::EventHeader>       evtHdr(new mu2e::EventHeader);
  std::unique_ptr<mu2e::EventWindowMarker> ewm((ewm_) ? new mu2e::EventWindowMarker : nullptr);
  art::Handle<artdaq::Fragments>           cfoFragmentHandle;

  if(!event.getByLabel(cfoFragmentTag_, cfoFragmentHandle)) {
    event.put(std::move(evtHdr));
    if(ewm_) event.put(std::move(ewm));
    TLOG(TLVL_DEBUG + 1) << "No CFO fragments found";
    return;
  }

  const auto *fragments = cfoFragmentHandle.product();
  if (fragments->size()>0){
    const auto &frag = fragments->at(0);
    const mu2e::CFOEventFragment   cfoFrag(frag);
    const CFOLib::CFO_Event        cfo       = cfoFrag.getData();
    const CFOLib::CFO_EventRecord& cfoRecord = cfo.GetEventRecord();
    evtHdr->mode          = 0;//cfo.GetEventMode();
    evtHdr->ewt           = static_cast<long unsigned int>(cfo.GetEventWindowTag().GetEventWindowTag().to_ullong());
    evtHdr->flags         = cfo.GetEventMode().isOnSpillFlagSet();
    evtHdr->eventDuration = cfoRecord.event_duration;
    TLOG(TLVL_DEBUG + 20) << "mode = " << evtHdr->mode << " ewt  = "<< evtHdr->ewt << " flags = " << evtHdr->flags
                          << " onspill = " << evtHdr->isOnSpill() << " duration = " << evtHdr->eventDuration;
  } else {
    TLOG(TLVL_DEBUG + 3) << "No CFO fragments found in the event";
  }

  if(ewm_) {
    constexpr double tick = 5.; // clock ticks -> ns
    ewm->_spillType   = (evtHdr->isOnSpill()) ? mu2e::EventWindowMarker::SpillType::onspill : mu2e::EventWindowMarker::SpillType::offspill;
    ewm->_eventLength = tick*evtHdr->eventDuration;
  }

  event.put(std::move(evtHdr));
  if(ewm_) event.put(std::move(ewm));
} // produce()

// ======================================================================

DEFINE_ART_MODULE(art::EventHeaderFromCFOFragment)

// ======================================================================
