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
#include "messagefacility/MessageLogger/MessageLogger.h"

// artdaq
#include <artdaq-core-mu2e/Data/EventHeader.hh>
#include "artdaq-core-mu2e/Overlays/CFOEventFragment.hh"
#include "artdaq-core-mu2e/Overlays/CFO_Packets/CFO_Event.h"
#include "artdaq-core-mu2e/Overlays/CFO_Packets/CFO_EventRecord.h"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include <artdaq-core/Data/Fragment.hh>

// Offline
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"

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
    fhicl::Atom<bool>          isSim    {Name("isSim")    , Comment("Flag to indicate simulation processing"), false};
    fhicl::Atom<int>           diagLevel{Name("diagLevel"), Comment("diagnostic level"), 0};
    fhicl::Atom<int>           maxWarnings{Name("maxWarnings"), Comment("Maximum number of missing CFO data warnings to print"), 10};
  };

  // --- C'tor/d'tor:
  explicit EventHeaderFromCFOFragment(const art::EDProducer::Table<Config>& config);
  ~EventHeaderFromCFOFragment() override {}

  void beginRun(art::Run&) override;
  void endJob() override;

  // --- Production:
  void produce(Event&) override;

private:
  art::InputTag cfoFragmentTag_;
  bool          ewm_;
  bool          isSim_;
  int           diagLevel_;
  int           maxWarnings_;
  // Counters for events without CFO data, reported at the end of the job
  unsigned long nEvents_      = 0;
  unsigned long nNoCFOHandle_ = 0;
  unsigned long nNoCFOFrags_  = 0;

  void warnMissing(const Event& event, unsigned long count, const std::string& msg) const;
};

// ======================================================================

void art::EventHeaderFromCFOFragment::beginRun(art::Run& Run) {}

// ----------------------------------------------------------------------

void art::EventHeaderFromCFOFragment::endJob() {
  if(nNoCFOHandle_ > 0 || nNoCFOFrags_ > 0) {
    mf::LogWarning("EventHeaderFromCFOFragment")
      << "Default EventHeader produced without CFO data in " << nNoCFOHandle_ + nNoCFOFrags_
      << " of " << nEvents_ << " events (" << nNoCFOHandle_ << " with no CFO fragment collection, "
      << nNoCFOFrags_ << " with an empty CFO fragment collection)";
  } else if(diagLevel_ > 0) {
    mf::LogInfo("EventHeaderFromCFOFragment") << "Found CFO data in all " << nEvents_ << " events";
  }
}

// ----------------------------------------------------------------------

void art::EventHeaderFromCFOFragment::warnMissing(const Event& event, unsigned long count,
                                                  const std::string& msg) const {
  if(maxWarnings_ < 0 || count <= static_cast<unsigned long>(maxWarnings_)) {
    mf::LogWarning("EventHeaderFromCFOFragment")
      << "Event " << event.id() << ": " << msg << " (" << cfoFragmentTag_.encode()
      << "), producing a default EventHeader"
      << ((count == static_cast<unsigned long>(maxWarnings_)) ? ". Further warnings suppressed" : "");
  }
}

art::EventHeaderFromCFOFragment::EventHeaderFromCFOFragment(
    const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    cfoFragmentTag_(config().cfoTag()),
    ewm_(config().ewm()),
    isSim_(config().isSim()),
    diagLevel_(config().diagLevel()),
    maxWarnings_(config().maxWarnings()) {
  produces<mu2e::EventHeader>();
  if(ewm_) {
    TLOG(TLVL_DEBUG + 2) << "Producing EventWindowMarker";
    produces<mu2e::EventWindowMarker>();
    produces<mu2e::ProtonBunchTime>();
  }
  if(isSim_) TLOG(TLVL_DEBUG + 2) << "Assuming simulation inputs!";
}

// ----------------------------------------------------------------------

void art::EventHeaderFromCFOFragment::produce(Event& event) {
  ++nEvents_;

  // Collection of CaloHits for the event
  std::unique_ptr<mu2e::EventHeader>       evtHdr(new mu2e::EventHeader);
  std::unique_ptr<mu2e::EventWindowMarker> ewm((ewm_) ? new mu2e::EventWindowMarker : nullptr);
  std::unique_ptr<mu2e::ProtonBunchTime>   pbt((ewm_) ? new mu2e::ProtonBunchTime   : nullptr);
  art::Handle<artdaq::Fragments>           cfoFragmentHandle;

  // Set values that are useful in simulations
  if(isSim_) {
    // Set the EWT to be the art event number for now
    evtHdr->ewt = static_cast<long unsigned int>(event.event());
  }

  // Check for the CFO fragment list
  if(!event.getByLabel(cfoFragmentTag_, cfoFragmentHandle)) {
    event.put(std::move(evtHdr));
    if(ewm_) {
      event.put(std::move(ewm));
      event.put(std::move(pbt));
    }
    TLOG(TLVL_DEBUG + 1) << "No CFO fragments found";
    warnMissing(event, ++nNoCFOHandle_, "No CFO fragment collection found");
    return;
  }

  // Process the first fragment found
  const auto *fragments = cfoFragmentHandle.product();
  if (fragments->size()>0){
    const auto &frag = fragments->at(0);
    const mu2e::CFOEventFragment   cfoFrag(frag);
    const CFOLib::CFO_Event        cfo       = cfoFrag.getData();
    const CFOLib::CFO_EventRecord& cfoRecord = cfo.GetEventRecord();
    evtHdr->mode          = cfoRecord.event_mode; //cfo.GetEventMode();
    evtHdr->ewt           = static_cast<long unsigned int>(cfo.GetEventWindowTag().GetEventWindowTag().to_ullong());
    evtHdr->flags         = cfo.GetEventMode().isOnSpillFlagSet();
    evtHdr->eventDuration = cfoRecord.event_duration;
    TLOG(TLVL_DEBUG + 20) << "mode = " << evtHdr->mode << " (" << cfoRecord.event_mode << ")"
                          << " ewt  = "<< evtHdr->ewt << " flags = " << int(evtHdr->flags)
                          << " onspill = " << evtHdr->isOnSpill() << " duration = " << evtHdr->eventDuration;
    if(diagLevel_ > 1) {
      mf::LogInfo("EventHeaderFromCFOFragment")
        << "Event " << event.id() << ": mode = " << evtHdr->mode << " ewt = " << evtHdr->ewt
        << " onspill = " << evtHdr->isOnSpill() << " duration = " << evtHdr->eventDuration;
    }
    if(ewm_) {
      constexpr double tick = 25.; // clock ticks -> ns
      ewm->_spillType   = (cfo.GetEventMode().isOnSpillFlagSet()) ? mu2e::EventWindowMarker::SpillType::onspill : mu2e::EventWindowMarker::SpillType::offspill;
      ewm->_eventLength = tick * evtHdr->eventDuration;
      pbt->pbtime_ = 0.f; // FIXME
      pbt->pbterr_ = 1.f;
    }

  } else {
    TLOG(TLVL_DEBUG + 3) << "No CFO fragments found in the event";
    warnMissing(event, ++nNoCFOFrags_, "CFO fragment collection is empty");
  }

  event.put(std::move(evtHdr));
  if(ewm_) {
    event.put(std::move(ewm));
    event.put(std::move(pbt));
  }
} // produce()

// ======================================================================

DEFINE_ART_MODULE(art::EventHeaderFromCFOFragment)

// ======================================================================
