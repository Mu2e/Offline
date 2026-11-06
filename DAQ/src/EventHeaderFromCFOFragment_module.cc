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
#include "Offline/DataProducts/inc/CFOEventTiming.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/DAQConditions/inc/EventTiming.hh"

// C++
#include <iostream>
#include <sstream>
#include <string>
#include <array>
#include <list>
#include <cmath>
#include <memory>
#include <optional>
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
  // On-spill events whose DR marker TDC is outside the valid range (current CFO firmware writes
  // placeholder values), so the nominal proton bunch time is used
  unsigned long nOnSpill_     = 0;
  unsigned long nBadTDC_      = 0;
  // Only needed for the EWM/PBT, so only set when producing them
  std::optional<mu2e::ProditionsHandle<mu2e::EventTiming>> eventTiming_h_;
  double clockTick_ = 0.; // 40 MHz clock tick in ns

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
  if(nBadTDC_ > 0) {
    mf::LogWarning("EventHeaderFromCFOFragment")
      << "DR marker TDC out of range in " << nBadTDC_ << " of " << nOnSpill_
      << " on-spill events: used the nominal proton bunch time instead of the TDC";
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
  if(hasInstance && cfoInstance_.find(':') != std::string::npos) {
    throw cet::exception("EventHeaderFromCFOFragment") << "cfoInstance tag ("
                                                       << cfoInstance_ << ") must not have a module:instance separator \":\"";
  }
  produces<mu2e::EventHeader>();
  if(ewm_) {
    TLOG(TLVL_DEBUG + 2) << "Producing EventWindowMarker";
    produces<mu2e::EventWindowMarker>();
    produces<mu2e::ProtonBunchTime>();
    eventTiming_h_.emplace();
    clockTick_ = mu2e::GlobalConstantsHandle<mu2e::PhysicsParams>()->getNominalDAQClockTick();
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
  // The default constructor leaves the check bit-fields uninitialized; zero them as the other constructors do
  evtHdr->rnr_check = evtHdr->ndtc_check = evtHdr->dtc_check = evtHdr->ewt_check = 0;

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
    // The on-spill flag is bit 0 of the mode4 byte (bits 32-39 of the 48-bit event mode word).
    // CFO_Event::GetEventMode() (artdaq-core-mu2e v9_09_00) converts the mode word with
    // DTC_EventMode(event_mode), which only fills mode0, so its isOnSpillFlagSet() is always false.
    // Read the flag from the record instead.
    const bool onSpill    = (cfoRecord.event_mode >> 32) & 0x1;
    evtHdr->mode          = cfoRecord.event_mode; //cfo.GetEventMode();
    evtHdr->ewt           = static_cast<long unsigned int>(cfo.GetEventWindowTag().GetEventWindowTag().to_ullong());
    evtHdr->flags         = onSpill ? mu2e::EventHeader::spillMask : 0;
    // The CFO records the event window length in clock ticks minus one, e.g. 3999 for a 4000-tick
    // (100 us) off-spill window
    evtHdr->eventDuration = mu2e::CFOEventTiming::windowTicks(cfoRecord.event_duration);
    evtHdr->rfmTDC_est      = cfoRecord.DR_marker_N_est;
    evtHdr->rfmTDC_measured = cfoRecord.DR_marker_N_meas;
    TLOG(TLVL_DEBUG + 20) << "mode = " << evtHdr->mode << " (" << cfoRecord.event_mode << ")"
                          << " ewt  = "<< evtHdr->ewt << " flags = " << int(evtHdr->flags)
                          << " onspill = " << evtHdr->isOnSpill() << " duration = " << evtHdr->eventDuration
                          << " DR TDC est = " << int(evtHdr->rfmTDC_est) << " meas = " << int(evtHdr->rfmTDC_measured);
    if(diagLevel_ > 1) {
      mf::LogInfo("EventHeaderFromCFOFragment")
        << "Event " << event.id() << ": mode = " << evtHdr->mode << " ewt = " << evtHdr->ewt
        << " onspill = " << evtHdr->isOnSpill() << " duration = " << evtHdr->eventDuration
        << " DR TDC est = " << int(evtHdr->rfmTDC_est) << " meas = " << int(evtHdr->rfmTDC_measured);
    }
    if(ewm_) {
      ewm->_spillType   = (onSpill) ? mu2e::EventWindowMarker::SpillType::onspill : mu2e::EventWindowMarker::SpillType::offspill;
      ewm->_eventLength = clockTick_ * evtHdr->eventDuration;
      if(ewm->_spillType == mu2e::EventWindowMarker::SpillType::onspill) ++nOnSpill_;
      // A valid TDC is 0 to 19 (one clock tick). The current CFO firmware writes placeholder TDC values
      // (e.g. 0xAA measured, 0xDD estimated), so the TDC is only used when it is in range.
      const bool validTDC = mu2e::CFOEventTiming::validTDC(evtHdr->rfmTDC_measured, clockTick_);
      if(ewm->_spillType == mu2e::EventWindowMarker::SpillType::onspill && validTDC) {
        // The DR marker N measurement is the time from a 40 MHz clock edge to the DR marker arrival
        // at the CFO (Mu2e-doc-55726). The event window starts a fixed number of clock ticks after that
        // edge, so the proton bunch time relative to the event window start is a constant plus the
        // measured TDC time. The constant matches the simulation convention (EventWindowMarkerProducer),
        // where the bunch arrives 0-25 ns before -timeFromProtonsToDRMarker. Take the TDC bin center.
        const mu2e::EventTiming& eventTiming = eventTiming_h_->get(event.id());
        pbt->pbtime_ = mu2e::CFOEventTiming::pbtimeFromTDC(evtHdr->rfmTDC_measured, eventTiming.timeFromProtonsToDRMarker(), clockTick_);
        pbt->pbterr_ = mu2e::CFOEventTiming::pbtimeFromTDCError();
      } else if(ewm->_spillType == mu2e::EventWindowMarker::SpillType::onspill) {
        // No valid TDC: use the nominal proton bunch time, the center of the 25 ns window it can be
        // in, as in simulation without a TDC (EventWindowMarkerProducer)
        const mu2e::EventTiming& eventTiming = eventTiming_h_->get(event.id());
        pbt->pbtime_ = mu2e::CFOEventTiming::nominalPBTime(eventTiming.timeFromProtonsToDRMarker(), clockTick_);
        pbt->pbterr_ = mu2e::CFOEventTiming::nominalPBTimeError(clockTick_);
        if(nBadTDC_++ == 0) {
          mf::LogWarning("EventHeaderFromCFOFragment")
            << "Event " << event.id() << ": DR marker TDC measurement " << int(evtHdr->rfmTDC_measured)
            << " is outside the valid range 0-" << mu2e::CFOEventTiming::maxTDC(clockTick_) << " (CFO record format version "
            << int(cfoRecord.record_format_version) << "), using the nominal proton bunch time."
            << " Further events are counted and reported at the end of the job";
        }
      } else { // no DR marker off spill
        pbt->pbtime_ = 0.f;
        pbt->pbterr_ = 0.f;
      }
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
