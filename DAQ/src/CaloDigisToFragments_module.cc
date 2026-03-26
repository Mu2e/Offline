// ======================================================================
//
// CaloDigisToFragments_module: Create DTCEVT artdaq::Fragment collection
// from CaloDigiCollection
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

#include "artdaq-core-mu2e/Overlays/Decoders/CalorimeterDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_Event.h"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_DataHeaderPacket.h"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_EventHeader.h"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_SubEventHeader.h"
#include "artdaq-core-mu2e/Overlays/DTCEventFragment.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include <artdaq-core/Data/Fragment.hh>

#include "Offline/CaloConditions/inc/CaloDAQMap.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/DataProducts/inc/CaloRawSiPMId.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {
struct CaloBlockData {
  uint8_t dtcID{0};
  uint8_t linkID{0};
  uint16_t transferByteCount{0};
  uint16_t packetCount{0};
  std::vector<mu2e::CalorimeterDataDecoder::CalorimeterHitDataPacket> hits;
  std::vector<std::vector<uint16_t>> waveforms;
};

constexpr uint8_t kFormatVersion = 1;
constexpr uint8_t kCaloRocsPerDtc = 6;
constexpr size_t kDtcEventHeaderBytes = 24;
constexpr size_t kDtcSubEventHeaderBytes = 48;
static_assert(sizeof(DTCLib::DTC_EventHeader) == kDtcEventHeaderBytes,
              "Unexpected DTC_EventHeader size");
static_assert(sizeof(DTCLib::DTC_SubEventHeader) == kDtcSubEventHeaderBytes,
              "Unexpected DTC_SubEventHeader size");

void writeBits(std::array<uint16_t, 6>& words, size_t startBit, size_t bitLength, uint32_t value) {
  // The decoder reads bits in a swapped-word layout, so we write into that same mapping.
  for (size_t relBit = 0; relBit < bitLength; ++relBit) {
    size_t bitIndex = startBit + relBit;
    size_t wordIndex = (bitIndex / 16) ^ 0x1;
    size_t bitOffset = 15 - (bitIndex % 16);
    uint16_t bit = (value >> (bitLength - relBit - 1)) & 0x1;
    if (bit != 0) {
      words[wordIndex] |= static_cast<uint16_t>(1u << bitOffset);
    }
  }
}

std::array<uint16_t, 6> encodeHitHeader(mu2e::CalorimeterDataDecoder::CalorimeterHitDataPacket const& hit) {
  std::array<uint16_t, 6> words{};
  writeBits(words, 0, 12, hit.Reserved1);
  writeBits(words, 12, 8, hit.BoardID);
  writeBits(words, 20, 3, hit.DetectorID);
  writeBits(words, 23, 5, hit.ChannelID);
  writeBits(words, 28, 16, hit.Time);
  writeBits(words, 44, 16, hit.InPayloadEventWindowTag);
  writeBits(words, 60, 12, hit.Baseline);
  writeBits(words, 72, 10, hit.IndexOfMaxDigitizerSample);
  writeBits(words, 82, 4, hit.ErrorFlags);
  writeBits(words, 86, 10, hit.NumberOfSamples);
  return words;
}

std::vector<uint8_t> encodeWaveform(std::vector<uint16_t> const& waveform) {
  if (waveform.empty()) {
    return {};
  }

  size_t sampleGroups = (waveform.size() + 3) / 4;
  size_t logicalWordCount = sampleGroups * 3;
  if (logicalWordCount % 2 != 0) {
    logicalWordCount++;
  }

  std::vector<uint16_t> logicalWords(logicalWordCount, 0);
  for (size_t i = 0; i < waveform.size(); i += 4) {
    uint16_t s0 = waveform[i] & 0x0FFF;
    uint16_t s1 = (i + 1 < waveform.size()) ? (waveform[i + 1] & 0x0FFF) : 0;
    uint16_t s2 = (i + 2 < waveform.size()) ? (waveform[i + 2] & 0x0FFF) : 0;
    uint16_t s3 = (i + 3 < waveform.size()) ? (waveform[i + 3] & 0x0FFF) : 0;

    size_t base = (i / 4) * 3;
    logicalWords[base] = static_cast<uint16_t>((s0 << 4) | (s1 >> 8));
    logicalWords[base + 1] = static_cast<uint16_t>((s1 << 8) | (s2 >> 4));
    logicalWords[base + 2] = static_cast<uint16_t>((s2 << 12) | s3);
  }

  // Hardware/decoder expects each 16-bit word pair swapped.
  std::vector<uint16_t> rawWords(logicalWords.size(), 0);
  for (size_t i = 0; i < logicalWords.size(); ++i) {
    rawWords[i ^ 0x1] = logicalWords[i];
  }

  // Each hit occupies 12-byte chunks: 12-byte header + packed 12-bit samples + chunk padding.
  size_t hitBytes = static_cast<size_t>(12 * std::ceil((12.0 + 1.5 * waveform.size()) / 12.0));
  size_t waveformBytes = hitBytes - 12;
  std::vector<uint8_t> rawBytes(waveformBytes, 0xFF);
  std::memcpy(rawBytes.data(), rawWords.data(), std::min(rawBytes.size(), rawWords.size() * sizeof(uint16_t)));
  return rawBytes;
}
} // namespace

namespace art {
class CaloDigisToFragments;
}

class art::CaloDigisToFragments : public EDProducer {
public:
  struct Config {
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level"), 0};
    fhicl::Atom<bool> useOfflineID{fhicl::Name("useOfflineID"),
                                   fhicl::Comment("Input CaloDigis use offline SiPM IDs"),
                                   true};
    fhicl::Atom<bool> skipFragmentOnSizeMismatch{
      fhicl::Name("skipFragmentOnSizeMismatch"),
      fhicl::Comment("Skip emitting fragment when packed bytes and event header size differ"),
      false};
    fhicl::Atom<art::InputTag> caloDigiTag{fhicl::Name("caloDigiTag"),
                                           fhicl::Comment("Input CaloDigiCollection"),
                                           art::InputTag("CaloDigisFromDTCEvents")};
  };

  explicit CaloDigisToFragments(const art::EDProducer::Table<Config>& config);
  virtual ~CaloDigisToFragments() {}

  virtual void produce(Event&);
  virtual void endJob();

private:
  mu2e::ProditionsHandle<mu2e::CaloDAQMap> calodaqconds_h_;

  int diagLevel_;
  bool useOfflineID_;
  bool skipFragmentOnSizeMismatch_;
  art::InputTag caloDigiTag_;

  long int total_events_;
  long int total_digis_;

  static size_t waveformMaximumIndex(const std::vector<uint16_t>& waveform);

  static void putBlockInEvent(DTCLib::DTC_Event& currentEvent, uint8_t dtcID,
                              DTCLib::DTC_DataBlock const& thisBlock);

  void buildDtcEventFromDigis(art::Event const& event, mu2e::CaloDigiCollection const& caloDigis,
                              DTCLib::DTC_Event& dtcEvent);
};

art::CaloDigisToFragments::CaloDigisToFragments(const art::EDProducer::Table<Config>& config)
    : art::EDProducer{config}
    , diagLevel_(config().diagLevel())
    , useOfflineID_(config().useOfflineID())
    , skipFragmentOnSizeMismatch_(config().skipFragmentOnSizeMismatch())
    , caloDigiTag_(config().caloDigiTag()) {
  produces<artdaq::Fragments>();
  total_events_ = 0;
  total_digis_ = 0;
}

size_t art::CaloDigisToFragments::waveformMaximumIndex(const std::vector<uint16_t>& waveform) {
  if (waveform.empty()) {
    return 0;
  }
  return std::distance(waveform.begin(), std::max_element(waveform.begin(), waveform.end()));
}

void art::CaloDigisToFragments::putBlockInEvent(DTCLib::DTC_Event& currentEvent, uint8_t dtcID,
                                                DTCLib::DTC_DataBlock const& thisBlock) {
  auto subEvt = currentEvent.GetSubEventByDTCID(dtcID, DTCLib::DTC_Subsystem_Calorimeter);
  if (subEvt == nullptr) {
    DTCLib::DTC_SubEvent newSubEvt;
    newSubEvt.SetEventWindowTag(currentEvent.GetEventWindowTag());
    newSubEvt.SetSourceDTC(dtcID, DTCLib::DTC_Subsystem_Calorimeter);
    newSubEvt.AddDataBlock(thisBlock);
    currentEvent.AddSubEvent(newSubEvt);
  } else {
    subEvt->AddDataBlock(thisBlock);
  }
}

void art::CaloDigisToFragments::buildDtcEventFromDigis(art::Event const& event,
                                                       mu2e::CaloDigiCollection const& caloDigis,
                                                       DTCLib::DTC_Event& dtcEvent) {
  auto const& calodaqconds = calodaqconds_h_.get(event.id());

  std::map<uint8_t, std::map<uint8_t, CaloBlockData>> byDtcAndLink;
  for (auto const& digi : caloDigis) {
    mu2e::CaloRawSiPMId rawId;
    uint16_t detectorID = 0;

    if (useOfflineID_) {
      mu2e::CaloSiPMId offId(static_cast<uint16_t>(digi.SiPMID()));
      rawId = calodaqconds.rawId(offId);
      detectorID = offId.detType();
    } else {
      rawId = mu2e::CaloRawSiPMId(static_cast<uint16_t>(digi.SiPMID()));
    }

    if (!rawId.isValid()) {
      continue;
    }

    auto const diracID = rawId.dirac();
    auto const dtcID = static_cast<uint8_t>(diracID / kCaloRocsPerDtc);
    auto const linkID = static_cast<uint8_t>(diracID % kCaloRocsPerDtc);
    auto const channelID = rawId.ROCchannel();
    auto& block = byDtcAndLink[dtcID][linkID];
    block.dtcID = dtcID;
    block.linkID = linkID;

    auto& hit = block.hits.emplace_back();
    hit.Reserved1 = 0xAAA;
    hit.BoardID = diracID;
    hit.DetectorID = detectorID;
    hit.ChannelID = channelID;
    hit.Time = static_cast<uint16_t>(digi.t0());
    hit.InPayloadEventWindowTag = 0;
    hit.Baseline = 0;
    hit.ErrorFlags = 0;

    std::vector<uint16_t> waveform;
    waveform.reserve(digi.waveform().size());
    for (auto const sample : digi.waveform()) {
      waveform.push_back(static_cast<uint16_t>(sample));
    }

    hit.NumberOfSamples = waveform.size();
    hit.IndexOfMaxDigitizerSample =
      waveform.empty() ? 0 : std::min<size_t>(digi.peakpos(), waveform.size() - 1);
    block.waveforms.emplace_back(std::move(waveform));
  }

  for (auto& [dtcID, linkMap] : byDtcAndLink) {
    for (uint8_t linkID = 0; linkID < kCaloRocsPerDtc; ++linkID) {
      auto& block = linkMap[linkID];
      block.dtcID = dtcID;
      block.linkID = linkID;

      if (block.hits.empty()) {
        block.transferByteCount = 16;
        block.packetCount = 0;
      } else {
        block.transferByteCount = 16;
        for (auto const& wf : block.waveforms) {
          block.transferByteCount +=
              static_cast<uint16_t>(12 * std::ceil((12.0 + 1.5 * wf.size()) / 12.0));
        }

        while (block.transferByteCount % 16 != 0) {
          block.transferByteCount++;
        }
        block.packetCount = (block.transferByteCount - 16) / 16;
      }

      DTCLib::DTC_DataBlock thisBlock(block.transferByteCount);
      if (thisBlock.blockPointer == nullptr) {
        continue;
      }
      std::fill(thisBlock.allocBytes->begin(), thisBlock.allocBytes->end(), 0xFF);

      auto const ts = DTCLib::DTC_EventWindowTag(static_cast<uint64_t>(event.event()));
      DTCLib::DTC_DataHeaderPacket hdr(static_cast<DTCLib::DTC_Link_ID>(block.linkID),
                                       block.packetCount,
                                       0,
                                       block.dtcID,
                                       DTCLib::DTC_Subsystem_Calorimeter,
                                       kFormatVersion,
                                       ts,
                                       0);
      auto hdrPkt = hdr.ConvertToDataPacket();
      hdrPkt.SetByte(5, static_cast<uint8_t>(((block.packetCount & 0x0700) >> 8) |
                     (DTCLib::DTC_Subsystem_Calorimeter << 5)));

      size_t pos = 0;
      std::memcpy(thisBlock.allocBytes->data() + pos, hdrPkt.GetData(), 16);
      pos += 16;

      if (!block.hits.empty()) {
        for (size_t i = 0; i < block.hits.size(); ++i) {
          auto rawHeader = encodeHitHeader(block.hits[i]);
          std::memcpy(thisBlock.allocBytes->data() + pos, rawHeader.data(), 12);
          pos += 12;

          auto rawWaveform = encodeWaveform(block.waveforms[i]);
          if (!rawWaveform.empty()) {
            std::memcpy(thisBlock.allocBytes->data() + pos, rawWaveform.data(), rawWaveform.size());
            pos += rawWaveform.size();
          }
        }
      }

      putBlockInEvent(dtcEvent, block.dtcID, thisBlock);
    }
  }
}

void art::CaloDigisToFragments::produce(Event& event) {
  total_events_++;

  std::unique_ptr<artdaq::Fragments> fragments(new artdaq::Fragments());
  auto caloDigiHandle = event.getHandle<mu2e::CaloDigiCollection>(caloDigiTag_);

  if (!caloDigiHandle.isValid()) {
    if (diagLevel_ > 0) {
      std::cout << "[CaloDigisToFragments::produce] No valid CaloDigiCollection found!" << std::endl;
    }
    event.put(std::move(fragments));
    return;
  }

  auto const& caloDigis = *caloDigiHandle;
  total_digis_ += caloDigis.size();

  DTCLib::DTC_Event dtcEvent;
  dtcEvent.SetEventWindowTag(DTCLib::DTC_EventWindowTag(static_cast<uint64_t>(event.event())));
  buildDtcEventFromDigis(event, caloDigis, dtcEvent);

  auto const eventBytes = dtcEvent.GetEventByteCount();
  if (eventBytes > 0 && dtcEvent.GetSubEventCount() > 0) {
    // Build the exact in-memory layout consumed by DTCEventFragment::SetupEvent.
    std::vector<uint8_t> packed;
    packed.reserve(eventBytes);

    auto const* evHdr = dtcEvent.GetHeader();
    packed.insert(packed.end(), reinterpret_cast<uint8_t const*>(evHdr),
            reinterpret_cast<uint8_t const*>(evHdr) + kDtcEventHeaderBytes);

    for (auto const& subEvt : dtcEvent.GetSubEvents()) {
      auto const* subHdr = subEvt.GetHeader();
      packed.insert(packed.end(), reinterpret_cast<uint8_t const*>(subHdr),
            reinterpret_cast<uint8_t const*>(subHdr) + kDtcSubEventHeaderBytes);

      for (auto const& block : subEvt.GetDataBlocks()) {
        auto const* blk = reinterpret_cast<uint8_t const*>(block.GetRawBufferPointer());
        packed.insert(packed.end(), blk, blk + block.byteSize);
      }
    }

    bool sizeMismatch = packed.size() != eventBytes;
    if (sizeMismatch) {
      if (skipFragmentOnSizeMismatch_) {
        if (diagLevel_ > 0) {
          std::cout << "[CaloDigisToFragments::produce] WARNING: packed size " << packed.size()
                    << " differs from header event size " << eventBytes
                    << " - skipping fragment" << std::endl;
        }
        event.put(std::move(fragments));
        return;
      }
      if (diagLevel_ > 0) {
        std::cout << "[CaloDigisToFragments::produce] WARNING: packed size " << packed.size()
                  << " differs from header event size " << eventBytes
                  << " - emitting fragment" << std::endl;
      }
    }

    {
      auto fragPtr = artdaq::Fragment::FragmentBytes(packed.size());
      fragPtr->setUserType(mu2e::FragmentType::DTCEVT);
      fragPtr->setSequenceID(static_cast<artdaq::Fragment::sequence_id_t>(event.event()));
      fragPtr->setFragmentID(0);
      fragPtr->setTimestamp(static_cast<artdaq::Fragment::timestamp_t>(event.event()));
      if (!packed.empty()) {
        std::memcpy(fragPtr->dataBeginBytes(), packed.data(), packed.size());
      }

      if (diagLevel_ > 0) {
        mu2e::DTCEventFragment testFragment(*fragPtr);
        auto testSubevents =
            testFragment.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Calorimeter);
        size_t decodedHits = 0;
        for (auto const& se : testSubevents) {
          mu2e::CalorimeterDataDecoder decoder(se);
          for (size_t i = 0; i < decoder.block_count(); ++i) {
            auto hitVec = decoder.GetCalorimeterHitData(i);
            if (hitVec != nullptr) {
              decodedHits += hitVec->size();
            }
          }
        }
        std::cout << "[CaloDigisToFragments::produce] Validation: " << testSubevents.size()
                  << " CALO subevent(s), " << decodedHits << " decodable hit(s) in output"
                  << std::endl;
      }

      fragments->emplace_back(std::move(*fragPtr));
    }
  }

  if (diagLevel_ > 0) {
    std::cout << "[CaloDigisToFragments::produce] Run " << event.run() << ", subrun "
              << event.subRun() << ", event " << event.event() << " has " << caloDigis.size()
              << " CaloDigis and created " << fragments->size() << " DTCEVT fragment(s)" << std::endl;
  }

  event.put(std::move(fragments));
}

void art::CaloDigisToFragments::endJob() {
  if (diagLevel_ > 0) {
    std::cout << "\n ----- [CaloDigisToFragments] Summary ----- " << std::endl;
    std::cout << "Total events: " << total_events_ << std::endl;
    std::cout << "Total CaloDigis processed: " << total_digis_ << std::endl;
  }
}

DEFINE_ART_MODULE(art::CaloDigisToFragments)
