// ======================================================================
//
// CFOFragmentFromSim_module: Create the CFO artdaq::Fragment collection
// from simulated event timing, so that the trigger and reconstruction
// can be run entirely from fragments on simulated events.
//
// The CFO event record (CFOLib::CFO_EventRecord) is filled from the
// simulated EventHeader made by EventWindowMarkerProducer: on-spill flag,
// event duration in 40 MHz clock ticks (recorded as the length minus one, as
// the CFO firmware does), and the DR marker TDC. Inputs made
// before EventWindowMarkerProducer produced an EventHeader can instead be
// used through the EventWindowMarker and ProtonBunchTimeMC, from which the
// same quantities are computed.
//
// EventHeaderFromCFOFragment decodes the output, e.g. with
//   cfoInstance : "CFO"   or   cfoTag : "CFOFragmentFromSim:CFO"
//
// The event window tag, fragment sequence ID and timestamp are set to the
// art event number, as in the *DigisToFragments modules, so the CFO and DTC
// fragments of an event carry the same event window tag.
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "artdaq-core-mu2e/Data/EventHeader.hh"
#include "artdaq-core-mu2e/Overlays/CFO_Packets/CFO_EventRecord.h"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include <artdaq-core/Data/Fragment.hh>

#include "Offline/DAQConditions/inc/EventTiming.hh"
#include "Offline/DataProducts/inc/CFOEventTiming.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <string>

namespace {
  // The CFO event record is seven 64-bit words (artdaq-core-mu2e CFO_EventRecord.h)
  constexpr size_t kCFORecordBytes = 56;
  static_assert(sizeof(CFOLib::CFO_EventRecord) == kCFORecordBytes, "Unexpected CFO_EventRecord size");

  // The event mode is 5 bytes (mode0-mode4) in the low 40 bits of the 48-bit event mode word,
  // with the on-spill flag in bit 0 of mode4 (DTC_EventMode::isOnSpillFlagSet, CFO_Event::SetEventMode)
  constexpr int      kModeFlagsShift = 32;
  constexpr uint64_t kOnSpillBit     = uint64_t(1) << kModeFlagsShift;
}

namespace mu2e {

  class CFOFragmentFromSim : public art::EDProducer {
  public:
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> eventHeaderTag{Name("eventHeaderTag"),
        Comment("Simulated EventHeader (from EventWindowMarkerProducer)"), art::InputTag("EWMProducer")};
      fhicl::Atom<bool> allowEWMFallback{Name("allowEWMFallback"),
        Comment("If no EventHeader is found, compute the CFO record from the EventWindowMarker and ProtonBunchTimeMC"), true};
      fhicl::Atom<art::InputTag> ewmTag{Name("ewmTag"),
        Comment("EventWindowMarker used when no EventHeader is found"), art::InputTag("EWMProducer")};
      fhicl::Atom<art::InputTag> pbtmcTag{Name("pbtmcTag"),
        Comment("ProtonBunchTimeMC used when no EventHeader is found"), art::InputTag("EWMProducer")};
      fhicl::Atom<std::string> instance{Name("instance"),
        Comment("Product instance name of the output fragment collection (the online CFO collection is daq:CFO)"), "CFO"};
      // Defaults match the CFO fragments in data (run 125472, November 2026)
      fhicl::Atom<unsigned> fragmentID{Name("fragmentID"), Comment("artdaq fragment ID of the CFO fragment"), 0};
      fhicl::Atom<unsigned> recordFormatVersion{Name("recordFormatVersion"),
        Comment("CFO event record format version written in the record"), 0x10};
      fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("diagnostic level"), 0};
    };

    explicit CFOFragmentFromSim(const art::EDProducer::Table<Config>& config);

    void produce(art::Event& event) override;
    void endJob() override;

  private:
    struct SimTiming {
      bool     onSpill  = false;
      uint32_t mode     = 0; // low 32 bits of the event mode (EventHeader::mode)
      uint16_t duration = 0; // 40 MHz clock ticks
      uint8_t  tdcEst   = 0;
      uint8_t  tdcMeas  = 0;
    };

    bool timingFromHeader(art::Event const& event, SimTiming& timing) const;
    void timingFromEWM(art::Event const& event, SimTiming& timing);

    art::InputTag eventHeaderTag_;
    bool          allowEWMFallback_;
    art::InputTag ewmTag_;
    art::InputTag pbtmcTag_;
    std::string   instance_;
    unsigned      fragmentID_;
    unsigned      recordFormatVersion_;
    int           diagLevel_;

    double period_; // 40 MHz clock tick in ns

    ProditionsHandle<EventTiming> eventTiming_h_;

    unsigned long nEvents_     = 0;
    unsigned long nFromHeader_ = 0;
    unsigned long nFromEWM_    = 0;
    unsigned long nOnSpill_    = 0;
  };

  // ----------------------------------------------------------------------

  CFOFragmentFromSim::CFOFragmentFromSim(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    eventHeaderTag_(config().eventHeaderTag()),
    allowEWMFallback_(config().allowEWMFallback()),
    ewmTag_(config().ewmTag()),
    pbtmcTag_(config().pbtmcTag()),
    instance_(config().instance()),
    fragmentID_(config().fragmentID()),
    recordFormatVersion_(config().recordFormatVersion()),
    diagLevel_(config().diagLevel()),
    period_(GlobalConstantsHandle<PhysicsParams>()->getNominalDAQClockTick())
  {
    if(recordFormatVersion_ > 0xFF) {
      throw cet::exception("CFOFragmentFromSim") << "recordFormatVersion " << recordFormatVersion_
                                                 << " does not fit in 8 bits";
    }
    if(fragmentID_ > std::numeric_limits<artdaq::Fragment::fragment_id_t>::max()) {
      throw cet::exception("CFOFragmentFromSim") << "fragmentID " << fragmentID_ << " is out of range";
    }
    produces<artdaq::Fragments>(instance_);
  }

  // ----------------------------------------------------------------------

  // Copy the timing from the simulated EventHeader
  bool CFOFragmentFromSim::timingFromHeader(art::Event const& event, SimTiming& timing) const {
    auto hdrH = event.getHandle<mu2e::EventHeader>(eventHeaderTag_);
    if(!hdrH.isValid()) return false;
    const auto& hdr = *hdrH;
    timing.onSpill  = hdr.isOnSpill();
    timing.mode     = hdr.mode;
    timing.duration = hdr.eventDuration;
    timing.tdcEst   = hdr.rfmTDC_est;
    timing.tdcMeas  = hdr.rfmTDC_measured;
    return true;
  }

  // ----------------------------------------------------------------------

  // Compute the timing for inputs without a simulated EventHeader, following EventWindowMarkerProducer:
  // the proton bunch arrives 0-25 ns (one clock tick) before -timeFromProtonsToDRMarker, and the DR
  // marker TDC measures where in that tick (CFOEventTiming)
  void CFOFragmentFromSim::timingFromEWM(art::Event const& event, SimTiming& timing) {
    auto ewmH   = event.getHandle<mu2e::EventWindowMarker>(ewmTag_);
    auto pbtmcH = event.getHandle<mu2e::ProtonBunchTimeMC>(pbtmcTag_);
    if(!ewmH.isValid() || !pbtmcH.isValid()) {
      throw cet::exception("CFOFragmentFromSim")
        << "Event " << event.id() << ": no EventHeader (" << eventHeaderTag_.encode() << ") and no "
        << (ewmH.isValid() ? "" : "EventWindowMarker (" + ewmTag_.encode() + ") ")
        << (pbtmcH.isValid() ? "" : "ProtonBunchTimeMC (" + pbtmcTag_.encode() + ") ")
        << "to compute the CFO record from";
    }
    const auto& ewm = *ewmH;
    timing.onSpill  = ewm.spillType() == EventWindowMarker::SpillType::onspill;
    timing.mode     = 0;
    timing.duration = static_cast<uint16_t>(std::lround(ewm.eventLength() / period_));
    if(timing.onSpill) {
      const EventTiming& eventTiming = eventTiming_h_.get(event.id());
      const double phaseShift = pbtmcH->pbtime_ + eventTiming.timeFromProtonsToDRMarker(); // in [-period, 0]
      const int tdc = CFOEventTiming::tdcFromPhase(phaseShift, period_);
      timing.tdcMeas = timing.tdcEst = static_cast<uint8_t>(tdc); // no prediction error simulated
    } else { // no DR marker off spill
      timing.tdcMeas = timing.tdcEst = 0;
    }
  }

  // ----------------------------------------------------------------------

  void CFOFragmentFromSim::produce(art::Event& event) {
    ++nEvents_;

    SimTiming timing;
    if(timingFromHeader(event, timing)) {
      ++nFromHeader_;
    } else if(allowEWMFallback_) {
      timingFromEWM(event, timing);
      ++nFromEWM_;
    } else {
      throw cet::exception("CFOFragmentFromSim")
        << "Event " << event.id() << ": no EventHeader found with tag " << eventHeaderTag_.encode()
        << " and allowEWMFallback is false";
    }
    if(timing.onSpill) ++nOnSpill_;

    // Fill the CFO event record. The 32-bit EventHeader::mode only holds the low 4 mode bytes,
    // so the flags byte (mode4) is set from the spill status
    CFOLib::CFO_EventRecord record;
    record.record_format_version = recordFormatVersion_;
    record.event_tag             = static_cast<uint64_t>(event.event());
    record.event_mode            = static_cast<uint64_t>(timing.mode) | (timing.onSpill ? kOnSpillBit : 0);
    // The CFO records the event window length minus one tick
    record.event_duration        = CFOEventTiming::recordedDuration(timing.duration);
    record.DR_marker_N_est       = timing.tdcEst;
    record.DR_marker_N_meas      = timing.tdcMeas;
    // Not simulated, left at 0: Linux time, DR marker N+1 estimate and measurement, TDC of marker N from spill start

    auto frag = artdaq::Fragment::FragmentBytes(kCFORecordBytes);
    frag->setUserType(mu2e::FragmentType::CFO);
    frag->setSequenceID(static_cast<artdaq::Fragment::sequence_id_t>(event.event()));
    frag->setFragmentID(static_cast<artdaq::Fragment::fragment_id_t>(fragmentID_));
    frag->setTimestamp(static_cast<artdaq::Fragment::timestamp_t>(event.event()));
    std::memcpy(frag->dataBeginBytes(), &record, kCFORecordBytes);

    if(diagLevel_ > 1) {
      mf::LogInfo("CFOFragmentFromSim")
        << "Event " << event.id() << ": ewt = " << record.event_tag << " onspill = " << timing.onSpill
        << " mode = 0x" << std::hex << record.event_mode << std::dec << " duration = " << timing.duration
        << " (recorded " << record.event_duration << ")"
        << " DR TDC est = " << int(timing.tdcEst) << " meas = " << int(timing.tdcMeas);
    }

    auto fragments = std::make_unique<artdaq::Fragments>();
    fragments->emplace_back(std::move(*frag));
    event.put(std::move(fragments), instance_);
  }

  // ----------------------------------------------------------------------

  void CFOFragmentFromSim::endJob() {
    if(diagLevel_ > 0 || (nFromEWM_ > 0 && nFromHeader_ > 0)) {
      mf::LogInfo("CFOFragmentFromSim")
        << "Made CFO fragments for " << nEvents_ << " events (" << nOnSpill_ << " on spill): "
        << nFromHeader_ << " from the EventHeader (" << eventHeaderTag_.encode() << "), "
        << nFromEWM_ << " from the EventWindowMarker (" << ewmTag_.encode() << ") and ProtonBunchTimeMC ("
        << pbtmcTag_.encode() << ")";
    }
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CFOFragmentFromSim)
