// ======================================================================
//
// StrawDigisToFragments_module: Create DTCEVT artdaq::Fragment collection
// from StrawDigiCollection (+ StrawDigiADCWaveformCollection)
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

#include "artdaq-core-mu2e/Overlays/Decoders/TrackerDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_Event.h"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_DataHeaderPacket.h"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_EventHeader.h"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_SubEventHeader.h"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include <artdaq-core/Data/Fragment.hh>

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/TrackerConditions/inc/TrackerPanelMap.hh"
#include "Offline/DataProducts/inc/StrawEnd.hh"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {
using TrackerDataPacket = mu2e::TrackerDataDecoder::TrackerDataPacket;
using TrackerADCPacket = mu2e::TrackerDataDecoder::TrackerADCPacket;

struct TrackerHitPayload {
  TrackerDataPacket mainPacket;
  std::vector<TrackerADCPacket> adcPackets;
  TrackerHitPayload() {
    mainPacket.StrawIndex = 0;
    mainPacket.TDC0A = 0;
    mainPacket.TDC0B = 0;
    mainPacket.TOT0 = 0;
    mainPacket.EWMCounter = 0;
    mainPacket.TDC1A = 0;
    mainPacket.TDC1B = 0;
    mainPacket.TOT1 = 0;
    mainPacket.ErrorFlags = 0;
    mainPacket.NumADCPackets = 0;
    mainPacket.PMP = 0;
    mainPacket.ADC00 = 0;
    mainPacket.ADC01A = 0;
    mainPacket.ADC01B = 0;
    mainPacket.ADC02 = 0;
    mainPacket.unused1 = 0;
  }
};

struct TrackerBlockData {
  uint8_t dtcID{0};
  uint8_t linkID{0};
  uint16_t transferByteCount{0};
  uint16_t packetCount{0};
  std::vector<TrackerHitPayload> hits;
};

constexpr uint8_t kFormatVersion = 1;
constexpr uint8_t kTrackerRocsPerDtc = 6;
constexpr size_t kDtcEventHeaderBytes = 24;
constexpr size_t kDtcSubEventHeaderBytes = 48;
static_assert(sizeof(DTCLib::DTC_EventHeader) == kDtcEventHeaderBytes,
              "Unexpected DTC_EventHeader size");
static_assert(sizeof(DTCLib::DTC_SubEventHeader) == kDtcSubEventHeaderBytes,
              "Unexpected DTC_SubEventHeader size");
static_assert(sizeof(TrackerDataPacket) == 16, "Unexpected TrackerDataPacket size");
static_assert(sizeof(TrackerADCPacket) == 16, "Unexpected TrackerADCPacket size");

void printTrackerPanelMap(mu2e::TrackerPanelMap const& panelMap) {
  size_t nrows = 0;
  std::cout << "[StrawDigisToFragments] TrackerPanelMap entries (online->offline):" << std::endl;
  for (uint32_t dtc = 0; dtc < mu2e::TrackerPanelMap::kMaxPlanes; ++dtc) {
    for (uint32_t link = 0; link < mu2e::StrawId::_npanels; ++link) {
      auto const* row = panelMap.panel_map_by_online_ind(dtc, link);
      if (row == nullptr) {
        continue;
      }
      ++nrows;
      std::cout << "  DTC=" << row->dtc() << " link=" << row->link()
                << " -> plane=" << row->uniquePlane()
                << " panel=" << row->panel()
                << " mnid=" << row->mnid()
                << " ppid=" << row->ppid()
                << " zface=" << row->zface() << std::endl;
    }
  }
  std::cout << "[StrawDigisToFragments] TrackerPanelMap rows printed: " << nrows << std::endl;
}

size_t waveformPacketCount(mu2e::StrawDigiADCWaveform const& waveform) {
  auto const& samples = waveform.samples();
  if (samples.size() < 3) {
    return 0;
  }
  // Match ArtBinaryPacketsFromDigis behavior: one tracker header packet carries
  // the first 3 ADC samples, then 12 samples per TrackerADCPacket.
  return static_cast<size_t>((samples.size() - 3) / 12);
}

TrackerHitPayload encodeTrackerHit(mu2e::StrawDigi const& digi,
                                   mu2e::StrawDigiADCWaveform const& waveform,
                                   uint16_t eventWindowCounter,
                                   uint16_t strawIndexWord) {
  TrackerHitPayload payload;
  payload.mainPacket.StrawIndex = strawIndexWord;
  payload.mainPacket.SetTDC0(digi.TDC(mu2e::StrawEnd::cal));
  payload.mainPacket.SetTDC1(digi.TDC(mu2e::StrawEnd::hv));
  payload.mainPacket.TOT0 = digi.TOT(mu2e::StrawEnd::cal);
  payload.mainPacket.TOT1 = digi.TOT(mu2e::StrawEnd::hv);
  payload.mainPacket.EWMCounter = eventWindowCounter & 0xF;
  payload.mainPacket.PMP = digi.PMP();
  payload.mainPacket.ErrorFlags = 0;
  payload.mainPacket.unused1 = 0;

  auto const& samples = waveform.samples();
  for (size_t i = 0; i < 3; ++i) {
    payload.mainPacket.SetWaveform(i, (i < samples.size() ? samples[i] : 0));
  }

  size_t const numADCPackets = waveformPacketCount(waveform);
  payload.mainPacket.NumADCPackets = numADCPackets;
  payload.adcPackets.reserve(numADCPackets);
  for (size_t i = 0; i < numADCPackets; ++i) {
    TrackerADCPacket packet;
    for (size_t j = 0; j < 12; ++j) {
      size_t const sampleIndex = 3 + i * 12 + j;
      packet.SetWaveform(j, (sampleIndex < samples.size() ? samples[sampleIndex] : 0));
    }
    payload.adcPackets.push_back(packet);
  }

  return payload;
}

} // namespace

namespace mu2e {
class StrawDigisToFragments;
}

class mu2e::StrawDigisToFragments : public art::EDProducer {
public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("diagnostic level"), 0};
    fhicl::Atom<bool> fallbackToOfflineWhenMapMissing{
      Name("fallbackToOfflineWhenMapMissing"),
      Comment("When TrackerPanelMap lookup fails, use offline plane/panel as DTC/link"),
      true};
    fhicl::Atom<bool> forceOfflineAddressing{
      Name("forceOfflineAddressing"),
      Comment("Ignore TrackerPanelMap and encode dtc/link/strawindex directly from offline StrawId"),
      false};
    fhicl::Atom<art::InputTag> strawDigiTag{Name("strawDigiTag"),
                                            Comment("Input StrawDigiCollection"),
                                            art::InputTag("makeSD")};
    fhicl::Atom<art::InputTag> strawDigiADCTag{Name("strawDigiADCTag"),
                                               Comment("Input StrawDigiADCWaveformCollection"),
                                               art::InputTag("makeSD")};
  };

  explicit StrawDigisToFragments(const art::EDProducer::Table<Config>& config);
  virtual ~StrawDigisToFragments() {}

  virtual void produce(art::Event&);
  virtual void endJob();

private:
  int diagLevel_;
  bool fallbackToOfflineWhenMapMissing_;
  bool forceOfflineAddressing_;
  art::InputTag strawDigiTag_;
  art::InputTag strawDigiADCTag_;
  bool printedTrackerPanelMap_;

  mu2e::ProditionsHandle<mu2e::TrackerPanelMap> trackerPanelMap_h_;

  long int total_events_;
  long int total_digis_;

  static void putBlockInEvent(DTCLib::DTC_Event& currentEvent, uint8_t dtcID,
                              DTCLib::DTC_DataBlock const& thisBlock);

  void buildDtcEventFromDigis(art::Event const& event,
                              mu2e::StrawDigiCollection const& strawDigis,
                              mu2e::StrawDigiADCWaveformCollection const& strawADCs,
                              DTCLib::DTC_Event& dtcEvent);
};

mu2e::StrawDigisToFragments::StrawDigisToFragments(const art::EDProducer::Table<Config>& config)
    : art::EDProducer{config}
    , diagLevel_(config().diagLevel())
    , fallbackToOfflineWhenMapMissing_(config().fallbackToOfflineWhenMapMissing())
    , forceOfflineAddressing_(config().forceOfflineAddressing())
    , strawDigiTag_(config().strawDigiTag())
    , strawDigiADCTag_(config().strawDigiADCTag())
    , printedTrackerPanelMap_(false) {
  produces<artdaq::Fragments>();
  total_events_ = 0;
  total_digis_ = 0;
}

void mu2e::StrawDigisToFragments::putBlockInEvent(DTCLib::DTC_Event& currentEvent, uint8_t dtcID,
                                                 DTCLib::DTC_DataBlock const& thisBlock) {
  auto subEvt = currentEvent.GetSubEventByDTCID(dtcID, DTCLib::DTC_Subsystem_Tracker);
  if (subEvt == nullptr) {
    DTCLib::DTC_SubEvent newSubEvt;
    newSubEvt.SetEventWindowTag(currentEvent.GetEventWindowTag());
    newSubEvt.SetSourceDTC(dtcID, DTCLib::DTC_Subsystem_Tracker);
    newSubEvt.AddDataBlock(thisBlock);
    currentEvent.AddSubEvent(newSubEvt);
  } else {
    subEvt->AddDataBlock(thisBlock);
  }
}

void mu2e::StrawDigisToFragments::buildDtcEventFromDigis(
    art::Event const& event, mu2e::StrawDigiCollection const& strawDigis,
    mu2e::StrawDigiADCWaveformCollection const& strawADCs, DTCLib::DTC_Event& dtcEvent) {
  auto const& trackerPanelMap = trackerPanelMap_h_.get(event.id());
  if (diagLevel_ > 1 && !printedTrackerPanelMap_) {
    printTrackerPanelMap(trackerPanelMap);
    printedTrackerPanelMap_ = true;
  }
  std::map<uint8_t, std::map<uint8_t, TrackerBlockData>> byDtcAndLink;

  size_t const nHits = std::min(strawDigis.size(), strawADCs.size());
  if (strawDigis.size() != strawADCs.size() && diagLevel_ > 0) {
    std::cout << "[StrawDigisToFragments::buildDtcEventFromDigis] WARNING: StrawDigi count "
              << strawDigis.size() << " differs from StrawDigiADCWaveform count "
              << strawADCs.size() << "; using first " << nHits << " entries" << std::endl;
  }

  for (size_t i = 0; i < nHits; ++i) {
    auto const& digi = strawDigis[i];
    auto const& waveform = strawADCs[i];

    auto const plane = digi.strawId().getPlane();
    auto const panel = digi.strawId().getPanel();

    if (panel >= kTrackerRocsPerDtc) {
      if (diagLevel_ > 0) {
        std::cout << "[StrawDigisToFragments::buildDtcEventFromDigis] WARNING: panel " << panel
                  << " outside expected ROC link range [0," << static_cast<int>(kTrackerRocsPerDtc)
                  << ") - skipping digi" << std::endl;
      }
      continue;
    }

    auto const* tpm = forceOfflineAddressing_ ? nullptr : trackerPanelMap.panel_map_by_offline_ind(plane, panel);

    uint8_t dtcID = 0;
    uint8_t linkID = 0;
    uint16_t strawIndexWord = 0;

    if (forceOfflineAddressing_) {
      dtcID = static_cast<uint8_t>(plane);
      linkID = static_cast<uint8_t>(panel);
      strawIndexWord = digi.strawId().asUint16();
    } else if (tpm != nullptr) {
      dtcID = static_cast<uint8_t>(tpm->dtc());
      linkID = static_cast<uint8_t>(tpm->link());
      auto const straw = static_cast<uint16_t>(digi.strawId().straw() & mu2e::StrawId::_strawmsk);
      auto const mnidWord = static_cast<uint16_t>(tpm->mnid() << mu2e::StrawId::_panelsft);
      strawIndexWord = static_cast<uint16_t>(mnidWord | straw);
    } else {
      if (!fallbackToOfflineWhenMapMissing_) {
        if (diagLevel_ > 0) {
          std::cout << "[StrawDigisToFragments::buildDtcEventFromDigis] WARNING: no TrackerPanelMap"
                    << " entry for offline plane=" << plane << " panel=" << panel
                    << " - skipping digi" << std::endl;
        }
        continue;
      }
      dtcID = static_cast<uint8_t>(plane);
      linkID = static_cast<uint8_t>(panel);
      strawIndexWord = digi.strawId().asUint16();
      if (diagLevel_ > 1) {
        std::cout << "[StrawDigisToFragments::buildDtcEventFromDigis] INFO: no TrackerPanelMap"
                  << " entry for offline plane=" << plane << " panel=" << panel
                  << " - using offline fallback DTC/link" << std::endl;
      }
    }

    if (linkID >= kTrackerRocsPerDtc) {
      if (diagLevel_ > 0) {
        std::cout << "[StrawDigisToFragments::buildDtcEventFromDigis] WARNING: mapped link "
                  << static_cast<int>(linkID) << " outside expected ROC link range [0,"
                  << static_cast<int>(kTrackerRocsPerDtc) << ") for plane=" << plane
                  << " panel=" << panel << " dtc=" << static_cast<int>(dtcID)
                  << " - skipping digi" << std::endl;
      }
      continue;
    }

    auto& block = byDtcAndLink[dtcID][linkID];
    block.dtcID = dtcID;
    block.linkID = linkID;

    uint16_t const ewmCounter = static_cast<uint16_t>(event.event() & 0xF);

    block.hits.emplace_back(encodeTrackerHit(digi, waveform, ewmCounter, strawIndexWord));
  }

  for (auto& [dtcID, linkMap] : byDtcAndLink) {
    // Emit all tracker ROC links in order so decoders expecting contiguous
    // ROC indices (0..5) can parse mixed-subdetector streams robustly.
    for (uint8_t linkID = 0; linkID < kTrackerRocsPerDtc; ++linkID) {
      auto& block = linkMap[linkID];
      block.dtcID = dtcID;
      block.linkID = linkID;

      // Calculate packet and byte counts based on actual hits
      // packetCount in DTC_DataHeaderPacket excludes the 16-byte data header
      // packet itself; it counts only hit payload packets.
      uint16_t numPackets = 0;
      for (auto const& hit : block.hits) {
        numPackets += 1;  // main packet per hit
        numPackets += static_cast<uint16_t>(hit.adcPackets.size());  // ADC packets per hit
      }

      block.packetCount = numPackets;
      block.transferByteCount = static_cast<uint16_t>((numPackets + 1) * 16);

      if (diagLevel_ > 1) {
        std::cout << "[buildDtcEventFromDigis] DTC " << static_cast<int>(dtcID)
                  << " Link " << static_cast<int>(linkID)
                  << " hits: " << block.hits.size()
                  << " numPackets: " << numPackets
                  << " transferByteCount: " << block.transferByteCount << std::endl;
      }

      DTCLib::DTC_DataBlock thisBlock(block.transferByteCount);
      if (thisBlock.blockPointer == nullptr) {
        if (diagLevel_ > 0) {
          std::cout << "[buildDtcEventFromDigis] Failed to allocate DTC_DataBlock with size "
                    << block.transferByteCount << std::endl;
        }
        continue;
      }
      std::fill(thisBlock.allocBytes->begin(), thisBlock.allocBytes->end(), 0xFF);

      auto const ts = DTCLib::DTC_EventWindowTag(static_cast<uint64_t>(event.event()));
      DTCLib::DTC_DataHeaderPacket hdr(static_cast<DTCLib::DTC_Link_ID>(block.linkID),
                                       block.packetCount,
                                       0,
                                       block.dtcID,
                                       DTCLib::DTC_Subsystem_Tracker,
                                       kFormatVersion,
                                       ts,
                                       0);
      auto hdrPkt = hdr.ConvertToDataPacket();

      size_t pos = 0;
      std::memcpy(thisBlock.allocBytes->data() + pos, hdrPkt.GetData(), 16);
      pos += 16;

      for (auto const& hit : block.hits) {
        std::memcpy(thisBlock.allocBytes->data() + pos, &hit.mainPacket, sizeof(TrackerDataPacket));
        pos += sizeof(TrackerDataPacket);

        if (!hit.adcPackets.empty()) {
          size_t const adcBytes = hit.adcPackets.size() * sizeof(TrackerADCPacket);
          std::memcpy(thisBlock.allocBytes->data() + pos, hit.adcPackets.data(), adcBytes);
          pos += adcBytes;
        }
      }

      if (diagLevel_ > 1) {
        std::cout << "  After filling: pos = " << pos << ", transferByteCount = " << block.transferByteCount << std::endl;
      }

      putBlockInEvent(dtcEvent, block.dtcID, thisBlock);
    }
  }
}

void mu2e::StrawDigisToFragments::produce(art::Event& event) {
  total_events_++;

  std::unique_ptr<artdaq::Fragments> fragments(new artdaq::Fragments());
  auto strawDigiHandle = event.getHandle<mu2e::StrawDigiCollection>(strawDigiTag_);
  auto strawADCsHandle = event.getHandle<mu2e::StrawDigiADCWaveformCollection>(strawDigiADCTag_);

  if (!strawDigiHandle.isValid() || !strawADCsHandle.isValid()) {
    if (diagLevel_ > 0) {
      std::cout << "[StrawDigisToFragments::produce] Missing input collections: StrawDigi valid="
                << strawDigiHandle.isValid() << ", StrawDigiADCWaveform valid="
                << strawADCsHandle.isValid() << std::endl;
    }
    event.put(std::move(fragments));
    return;
  }

  auto const& strawDigis = *strawDigiHandle;
  auto const& strawADCs = *strawADCsHandle;
  total_digis_ += strawDigis.size();

  DTCLib::DTC_Event dtcEvent;
  dtcEvent.SetEventWindowTag(DTCLib::DTC_EventWindowTag(static_cast<uint64_t>(event.event())));
  buildDtcEventFromDigis(event, strawDigis, strawADCs, dtcEvent);

  auto const eventBytes = dtcEvent.GetEventByteCount();
  if (diagLevel_ > 1) {
    std::cout << "[StrawDigisToFragments::produce] eventBytes from DTC_Event: " << eventBytes
              << ", subEventCount: " << dtcEvent.GetSubEventCount() << std::endl;
    for (size_t i = 0; i < dtcEvent.GetSubEvents().size(); ++i) {
      auto const& subEvt = dtcEvent.GetSubEvents()[i];
      std::cout << "  SubEvent " << i << " blockCount: " << subEvt.GetDataBlocks().size() << std::endl;
      for (size_t j = 0; j < subEvt.GetDataBlocks().size(); ++j) {
        auto const& block = subEvt.GetDataBlocks()[j];
        std::cout << "    Block " << j << " byteSize: " << block.byteSize << std::endl;
      }
    }
  }

  if (dtcEvent.GetSubEventCount() > 0) {
    // Emit one fragment per tracker subevent (per DTC). The downstream decoder
    // processes a single subevent per artdaq fragment payload.
    for (size_t seIdx = 0; seIdx < dtcEvent.GetSubEvents().size(); ++seIdx) {
      auto const& subEvt = dtcEvent.GetSubEvents()[seIdx];

      size_t calcEventBytes = kDtcEventHeaderBytes + kDtcSubEventHeaderBytes;
      for (auto const& block : subEvt.GetDataBlocks()) {
        calcEventBytes += block.byteSize;
      }

      std::vector<uint8_t> packed;
      packed.reserve(calcEventBytes);

      auto const* evHdr = dtcEvent.GetHeader();
      packed.insert(packed.end(), reinterpret_cast<uint8_t const*>(evHdr),
                    reinterpret_cast<uint8_t const*>(evHdr) + kDtcEventHeaderBytes);

      auto const* subHdr = subEvt.GetHeader();
      packed.insert(packed.end(), reinterpret_cast<uint8_t const*>(subHdr),
                    reinterpret_cast<uint8_t const*>(subHdr) + kDtcSubEventHeaderBytes);

      for (auto const& block : subEvt.GetDataBlocks()) {
        auto const* blk = reinterpret_cast<uint8_t const*>(block.GetRawBufferPointer());
        packed.insert(packed.end(), blk, blk + block.byteSize);
      }

      // Each output fragment carries exactly one subevent payload. Patch the
      // copied event header so generic DTC decoders see consistent metadata.
      auto* packedEventHdr = reinterpret_cast<DTCLib::DTC_EventHeader*>(packed.data());
      packedEventHdr->inclusive_event_byte_count = static_cast<uint32_t>(calcEventBytes);
      packedEventHdr->num_dtcs = 1;

      if (diagLevel_ > 1 && packed.size() != calcEventBytes) {
        std::cout << "[StrawDigisToFragments::produce] WARNING: subevent " << seIdx
                  << " packed size " << packed.size()
                  << " differs from recalculated " << calcEventBytes << std::endl;
      }

      auto fragPtr = artdaq::Fragment::FragmentBytes(packed.size());
      fragPtr->setUserType(mu2e::FragmentType::DTCEVT);
      fragPtr->setSequenceID(static_cast<artdaq::Fragment::sequence_id_t>(event.event()));
      fragPtr->setFragmentID(static_cast<artdaq::Fragment::fragment_id_t>(seIdx));
      fragPtr->setTimestamp(static_cast<artdaq::Fragment::timestamp_t>(event.event()));
      if (!packed.empty()) {
        std::memcpy(fragPtr->dataBeginBytes(), packed.data(), packed.size());
      }

      fragments->emplace_back(std::move(*fragPtr));
    }

    if (diagLevel_ > 0 && eventBytes != 0) {
      size_t aggregatePacked = 0;
      for (auto const& subEvt : dtcEvent.GetSubEvents()) {
        aggregatePacked += kDtcEventHeaderBytes + kDtcSubEventHeaderBytes;
        for (auto const& block : subEvt.GetDataBlocks()) aggregatePacked += block.byteSize;
      }
      if (aggregatePacked != eventBytes * dtcEvent.GetSubEventCount()) {
        std::cout << "[StrawDigisToFragments::produce] INFO: GetEventByteCount()=" << eventBytes
                  << " is unreliable for multi-subevent packing; emitted "
                  << fragments->size() << " per-subevent fragment(s)" << std::endl;
      }
    }
  }

  if (diagLevel_ > 0) {
    std::cout << "[StrawDigisToFragments::produce] Run " << event.run() << ", subrun "
              << event.subRun() << ", event " << event.event() << " has " << strawDigis.size()
              << " StrawDigis and created " << fragments->size() << " DTCEVT fragment(s)"
              << std::endl;
  }

  event.put(std::move(fragments));
}

void mu2e::StrawDigisToFragments::endJob() {
  if (diagLevel_ > 0) {
    std::cout << "\n ----- [StrawDigisToFragments] Summary ----- " << std::endl;
    std::cout << "Total events: " << total_events_ << std::endl;
    std::cout << "Total StrawDigis processed: " << total_digis_ << std::endl;
  }
}

DEFINE_ART_MODULE(mu2e::StrawDigisToFragments)
