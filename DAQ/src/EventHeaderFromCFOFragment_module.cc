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
#include "fhiclcpp/types/OptionalAtom.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "cetlib_except/exception.h"
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
#include <sstream>
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
    fhicl::OptionalAtom<art::InputTag> cfoTag     {Name("cfoTag")     , Comment("Input CFO fragment tag (set this or cfoInstance)")};
    fhicl::OptionalAtom<std::string>   cfoInstance{Name("cfoInstance"), Comment("Product instance name of the CFO fragment collection (e.g. CFO), from any module label or process (e.g. daq:CFO or aggregator:CFO); set this or cfoTag")};
    fhicl::Atom<bool>        ewm        {Name("createEWM")  , Comment("Produce an event window marker object")};
    fhicl::Atom<bool>        isSim      {Name("isSim")      , Comment("Simulation processing: no CFO data is expected, so missing CFO data is not a warning"), false};
    fhicl::Atom<int>         diagLevel  {Name("diagLevel")  , Comment("diagnostic level"), 0};
  };

  // --- C'tor/d'tor:
  explicit EventHeaderFromCFOFragment(const art::EDProducer::Table<Config>& config);
  ~EventHeaderFromCFOFragment() override {}

  void beginRun(art::Run&) override;
  void endJob() override;

  // --- Production:
  void produce(Event&) override;

private:
  art::InputTag cfoTag_;
  std::string cfoInstance_; // if set, select by product instance name instead of cfoTag_
  bool        ewm_;
  bool        isSim_;
  int         diagLevel_;
  // Counters for events without CFO data, reported at the end of the job
  unsigned long nEvents_      = 0;
  unsigned long nNoCFOHandle_ = 0;
  unsigned long nNoCFOFrags_  = 0;

  void warnMissing(const Event& event, const std::string& msg) const;
  std::string cfoSource() const { return cfoInstance_.empty() ? cfoTag_.encode() : "instance " + cfoInstance_; }
};

// ======================================================================

void art::EventHeaderFromCFOFragment::beginRun(art::Run& Run) {}

// ----------------------------------------------------------------------

void art::EventHeaderFromCFOFragment::endJob() {
  if(nNoCFOHandle_ > 0 || nNoCFOFrags_ > 0) {
    // Simulation has no CFO data, so the default header is the expected result there
    if(isSim_ && diagLevel_ <= 0) return;
    std::ostringstream msg;
    msg << "Default EventHeader produced without CFO data in " << nNoCFOHandle_ + nNoCFOFrags_
        << " of " << nEvents_ << " events (" << nNoCFOHandle_ << " with no CFO fragment collection, "
        << nNoCFOFrags_ << " with an empty CFO fragment collection)";
    if(isSim_) mf::LogInfo   ("EventHeaderFromCFOFragment") << msg.str();
    else       mf::LogWarning("EventHeaderFromCFOFragment") << msg.str();
  } else if(diagLevel_ > 0) {
    mf::LogInfo("EventHeaderFromCFOFragment") << "Found CFO data in all " << nEvents_ << " events";
  }
}

// ----------------------------------------------------------------------

// Repeated warnings are rate-limited by the message facility configuration for this category.
// In simulation no CFO data is expected, so this is only reported at diagLevel > 1.
void art::EventHeaderFromCFOFragment::warnMissing(const Event& event, const std::string& msg) const {
  if(isSim_) {
    if(diagLevel_ > 1) {
      mf::LogInfo("EventHeaderFromCFOFragment")
        << "Event " << event.id() << ": " << msg << " (" << cfoSource()
        << "), producing a default EventHeader";
    }
    return;
  }
  mf::LogWarning("EventHeaderFromCFOFragment")
    << "Event " << event.id() << ": " << msg << " (" << cfoSource()
    << "), producing a default EventHeader";
}

art::EventHeaderFromCFOFragment::EventHeaderFromCFOFragment(
    const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    ewm_(config().ewm()),
    isSim_(config().isSim()),
    diagLevel_(config().diagLevel()) {
  const bool hasTag      = config().cfoTag(cfoTag_);
  const bool hasInstance = config().cfoInstance(cfoInstance_);
  if(hasTag == hasInstance) {
    throw cet::exception("EventHeaderFromCFOFragment")
      << "Exactly one of cfoTag (e.g. daq:CFO) or cfoInstance (e.g. CFO) must be set";
  }
  if(hasInstance && cfoInstance_.empty()) {
    throw cet::exception("EventHeaderFromCFOFragment") << "cfoInstance must not be empty";
  }
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
  // Value-initialize so the defaults put without CFO data are well defined.
  // The EWT is left at 0 in simulation for now.
  std::unique_ptr<mu2e::EventHeader>       evtHdr(new mu2e::EventHeader);
  std::unique_ptr<mu2e::EventWindowMarker> ewm((ewm_) ? new mu2e::EventWindowMarker{} : nullptr);
  std::unique_ptr<mu2e::ProtonBunchTime>   pbt((ewm_) ? new mu2e::ProtonBunchTime{}   : nullptr);
  art::Handle<artdaq::Fragments>           cfoFragmentHandle;

  // Check for the CFO fragment list
  if(!cfoInstance_.empty()) { // select by product instance name only
    const auto handles = event.getMany<artdaq::Fragments>(art::ProductInstanceNameSelector{cfoInstance_});
    if(handles.size() > 1) {
      throw cet::exception("EventHeaderFromCFOFragment")
        << "Event " << event.id() << ": " << handles.size() << " fragment collections have product instance name '"
        << cfoInstance_ << "', expected at most one";
    }
    if(handles.size() == 1 && handles[0].isValid()) cfoFragmentHandle = handles[0];
  } else { // get this specific collection
    event.getByLabel(cfoTag_, cfoFragmentHandle);
  }

  if(!cfoFragmentHandle.isValid()) {
    event.put(std::move(evtHdr));
    if(ewm_) {
      event.put(std::move(ewm));
      event.put(std::move(pbt));
    }
    TLOG(TLVL_DEBUG + 1) << "No CFO fragments found";
    ++nNoCFOHandle_;
    warnMissing(event, "No CFO fragment collection found");
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
    ++nNoCFOFrags_;
    warnMissing(event, "CFO fragment collection is empty");
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
