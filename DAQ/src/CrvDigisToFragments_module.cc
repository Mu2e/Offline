// ======================================================================
//
// CrvDigisToFragments_module: Create DTCEVT artdaq::Fragment collection
// from CrvDigiCollection (+ CrvDigiADCWaveformCollection)
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

#include "artdaq-core-mu2e/Overlays/Decoders/CRVDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_Event.h"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_DataHeaderPacket.h"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_EventHeader.h"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_SubEventHeader.h"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include "artdaq-core-mu2e/Overlays/DTCEventFragment.hh"
#include <artdaq-core/Data/Fragment.hh>

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/CRVConditions/inc/CRVOrdinal.hh"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace mu2e
{

constexpr uint8_t kFormatVersion = 1;
constexpr uint8_t kBytesPerPacket = 16;
constexpr size_t kDtcEventHeaderBytes = 24;
constexpr size_t kDtcSubEventHeaderBytes = 48;
static_assert(sizeof(DTCLib::DTC_EventHeader) == kDtcEventHeaderBytes, "Unexpected DTC_EventHeader size");
static_assert(sizeof(DTCLib::DTC_SubEventHeader) == kDtcSubEventHeaderBytes, "Unexpected DTC_SubEventHeader size");

struct CrvBlockData
{
  CRVDataDecoder::CRVROCStatusPacketFEBII crvROCstatus;
  std::vector<CRVDataDecoder::CRVHitRawFEBII> crvHits;
};

CRVDataDecoder::CRVHitRawFEBII encodeCrvHit(CrvDigi const& digi)
{
  CRVDataDecoder::CRVHitRawFEBII hitPacket;
  hitPacket.hitInfo.portNumber  = digi.GetFEB();
  hitPacket.hitInfo.fpgaNumber  = digi.GetFEBchannel()/CRVId::nChanPerFPGA;
  hitPacket.hitInfo.fpgaChannel = digi.GetFEBchannel()%CRVId::nChanPerFPGA;
  hitPacket.hitInfo.hitTime     = digi.GetStartTDC()*2; //online: period of 6.25ns, offline: period of 12.5ns.

  const std::vector<int16_t> &adcs = digi.GetADCs();  //check for correct number of samples in calling function
  for(size_t adcBlock=0; adcBlock<CRVDataDecoder::nADCblocks; ++adcBlock)
  {
    hitPacket.adcBlocks[adcBlock].ADCsample0  = adcs.at(adcBlock*CRVDataDecoder::nADCsamplesPerBlock+0);
    uint16_t ADCsample1 = adcs.at(adcBlock*CRVDataDecoder::nADCsamplesPerBlock+1);
    hitPacket.adcBlocks[adcBlock].ADCsample1a = ADCsample1 & 0xf; //lower 4 bit
    hitPacket.adcBlocks[adcBlock].ADCsample1b = (ADCsample1 & 0xff0) >> 4; //upper 8 bit
    uint16_t ADCsample2 = adcs.at(adcBlock*CRVDataDecoder::nADCsamplesPerBlock+2);
    hitPacket.adcBlocks[adcBlock].ADCsample2a = ADCsample2 & 0xff; //lower 8 bit
    hitPacket.adcBlocks[adcBlock].ADCsample2b = (ADCsample2 & 0xf00) >> 8; //upper 4 bit
    hitPacket.adcBlocks[adcBlock].ADCsample3  = adcs.at(adcBlock*CRVDataDecoder::nADCsamplesPerBlock+3);
  }

  return hitPacket;
}

class CrvDigisToFragments : public art::EDProducer
{
  public:
  struct Config
  {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<int>           diagLevel{Name("diagLevel"), Comment("diagnostic level"), 0};
    fhicl::Atom<art::InputTag> crvDigiTag{Name("CrvDigiTag"), Comment("Input CrvDigiCollection"), art::InputTag("makeSD")};
    fhicl::Atom<int>           crvDtcIdStart{Name("CrvDtcIdStart"), Comment("First CRV DTC ID"), 0};  //FIXME: Which DTC IDs do the CRV use?
  };

  explicit CrvDigisToFragments(const art::EDProducer::Table<Config>& config);
  virtual ~CrvDigisToFragments() {}

  virtual void produce(art::Event&);
  virtual void endJob();

  private:
  int           _diagLevel;
  art::InputTag _crvDigiTag;
  int           _crvDtcIdStart;

  long int _totalEvents;
  long int _totalDigis;

  std::set<uint16_t> _ROCs;

  ProditionsHandle<CRVOrdinal> _crvChannelMap_h;

  const int _fragmentIdOffset = 0; // IDs start from here

  static void putBlockInEvent(DTCLib::DTC_Event& currentEvent, uint8_t dtcID,
                              DTCLib::DTC_DataBlock const& thisBlock);

  void buildDtcEventsFromDigis(art::Event const& event,
                               mu2e::CrvDigiCollection const& crvDigis,
                               std::map<uint8_t, DTCLib::DTC_Event>& dtcEvents);
};

CrvDigisToFragments::CrvDigisToFragments(const art::EDProducer::Table<Config>& config)
    : art::EDProducer{config}
    , _diagLevel(config().diagLevel())
    , _crvDigiTag(config().crvDigiTag())
    , _crvDtcIdStart(config().crvDtcIdStart())
{
  produces<artdaq::Fragments>();
  _totalEvents = 0;
  _totalDigis = 0;
}

void CrvDigisToFragments::putBlockInEvent(DTCLib::DTC_Event& currentEvent, uint8_t dtcID,
                                          DTCLib::DTC_DataBlock const& thisBlock)
{
  auto subEvt = currentEvent.GetSubEventByDTCID(dtcID, DTCLib::DTC_Subsystem_CRV);
  if (subEvt == nullptr)
  {
    DTCLib::DTC_SubEvent newSubEvt;
    newSubEvt.SetEventWindowTag(currentEvent.GetEventWindowTag());
    newSubEvt.SetSourceDTC(dtcID, DTCLib::DTC_Subsystem_CRV);
    newSubEvt.AddDataBlock(thisBlock);
    currentEvent.AddSubEvent(newSubEvt);
  }
  else
  {
    subEvt->AddDataBlock(thisBlock);
    currentEvent.UpdateHeader();  //update the byte count
  }
}

void CrvDigisToFragments::buildDtcEventsFromDigis(art::Event const& event, mu2e::CrvDigiCollection const& crvDigis,
                                                  std::map<uint8_t, DTCLib::DTC_Event>& dtcEvents)
{
  std::map<uint8_t, std::map<uint8_t, CrvBlockData>> byDtcAndLink;

  for(size_t i = 0; i < crvDigis.size(); ++i)
  {
    auto const& digi = crvDigis[i];

    if(digi.IsNZS()) continue; //FEBs can't handle NZS data, yet

    uint16_t ROC = digi.GetROC();
    if(ROC==0) //shouldn't happen
    {
      if(_diagLevel > 0)
      {
        std::cout << "[CrvDigisToFragments::buildDtcEventsFromDigis] WARNING: Invalid ROC - skipping digi" << std::endl;
      }
      continue;
    }

    if(digi.GetADCs().size()!=CRVDataDecoder::nADCsamples)
    {
      if(_diagLevel > 0)
      {
        std::cout << "[CrvDigisToFragments::buildDrcEventsFromDigis] WARNING: Encountered waveform with " << digi.GetADCs().size() << " samples. "
                  << "Can handle only " << CRVDataDecoder::nADCsamples << " samples. - skipping digi" << std::endl;
      }
      continue;
    }

    uint8_t dtcID  = _crvDtcIdStart + (ROC-1)/CRVId::nROCPerDTC;
    uint8_t linkID = (ROC-1)%CRVId::nROCPerDTC;

    auto& block = byDtcAndLink[dtcID][linkID];

    block.crvHits.emplace_back(encodeCrvHit(digi));
  }

  for(auto& [dtcID, linkMap] : byDtcAndLink)
  {
    auto& dtcEvent = dtcEvents[dtcID];
    dtcEvent.SetEventWindowTag(DTCLib::DTC_EventWindowTag(static_cast<uint64_t>(event.event())));

    for(uint8_t linkID = 0; linkID < CRVId::nROCPerDTC; ++linkID)
    {
      auto& block = linkMap[linkID];

      // Calculate packet and byte counts based on actual hits
      // packetCount in DTC_DataHeaderPacket excludes the 16-byte data header
      // packet itself; it counts only hit payload packets.
      size_t byteCount = block.crvHits.size()*sizeof(CRVDataDecoder::CRVHitRawFEBII) + sizeof(CRVDataDecoder::CRVROCStatusPacketFEBII);
      size_t packetCount = (byteCount+kBytesPerPacket-1) / kBytesPerPacket;   //integer ceiling division of byteCount/bytesPerPacket
      size_t transferByteCount = (packetCount + 1) * kBytesPerPacket;
      DTCLib::DTC_Subsystem subsystem = DTCLib::DTC_Subsystem_CRV;

      uint16_t ROC = (dtcID-_crvDtcIdStart)*CRVId::nROCPerDTC + linkID + 1;
      if(_ROCs.find(ROC)==_ROCs.end()) //not a used ROC
      {
        byteCount=0;
        packetCount=0;
        transferByteCount=kBytesPerPacket;
        subsystem=static_cast<DTCLib::DTC_Subsystem>(0);
      }

      if(byteCount%2!=0)
      {
        if(_diagLevel > 0)
        {
          std::cout << "[buildDtcEventsFromDigis] Byte count is not a multiple of 2" << std::endl;
        }
        continue;
      }

      block.crvROCstatus.ControllerEventWordCount=byteCount/2;
      block.crvROCstatus.EventWindowTag0=event.event() & 0xffff;
      block.crvROCstatus.EventWindowTag1=(event.event()>>16) & 0xffff;

      if(_diagLevel > 1)
      {
        std::cout << "[buildDtcEventsFromDigis] DTC " << static_cast<int>(dtcID)
                  << " Link " << static_cast<int>(linkID)
                  << " hits: " << block.crvHits.size()
                  << " numPackets: " << packetCount
                  << " transferByteCount: " << transferByteCount << std::endl;
      }

      DTCLib::DTC_DataBlock thisBlock(transferByteCount);
      if(thisBlock.blockPointer == nullptr)
      {
        if(_diagLevel > 0)
        {
          std::cout << "[buildDtcEventsFromDigis] Failed to allocate DTC_DataBlock with size "
                    << transferByteCount << std::endl;
        }
        continue;
      }
      std::fill(thisBlock.allocBytes->begin(), thisBlock.allocBytes->end(), 0xFF);

      auto const ts = DTCLib::DTC_EventWindowTag(static_cast<uint64_t>(event.event()));
      DTCLib::DTC_DataHeaderPacket hdr(static_cast<DTCLib::DTC_Link_ID>(linkID),
                                       packetCount,
                                       0,
                                       dtcID,
                                       subsystem,
                                       kFormatVersion,
                                       ts,
                                       0);

      auto hdrPkt = hdr.ConvertToDataPacket();  //this function has a bug. it does not transfer the subsystem ID to the DTC_DataPacket.
      hdrPkt.SetByte(5, static_cast<uint8_t>(((packetCount & 0x0700) >> 8) | (subsystem << 5)));  //this is a temporary bug fix.

      size_t pos = 0;
      std::memcpy(thisBlock.allocBytes->data() + pos, hdrPkt.GetData(), kBytesPerPacket);
      pos += kBytesPerPacket;

      if(packetCount>0) //no ROC status packet, if ROC is not used
      {
        std::memcpy(thisBlock.allocBytes->data() + pos, &block.crvROCstatus, sizeof(CRVDataDecoder::CRVROCStatusPacketFEBII));
        pos += sizeof(CRVDataDecoder::CRVROCStatusPacketFEBII);

        for(auto const& hit : block.crvHits)
        {
          std::memcpy(thisBlock.allocBytes->data() + pos, &hit, sizeof(CRVDataDecoder::CRVHitRawFEBII));
          pos += sizeof(CRVDataDecoder::CRVHitRawFEBII);
        }
      }

      if(_diagLevel > 1)
      {
        std::cout << "  After filling: pos = " << pos << ", transferByteCount = " << transferByteCount << std::endl;
      }

      putBlockInEvent(dtcEvent, dtcID, thisBlock);
    }
  }
}

void CrvDigisToFragments::produce(art::Event& event)
{
  if(_ROCs.empty())
  {
    auto const& crvChannelMap = _crvChannelMap_h.get(event.id());
    GeomHandle<CosmicRayShield> CRS;
    const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
    for(size_t channel=0; channel<counters.size()*CRVId::nChanPerBar; ++channel)
    {
      if(!crvChannelMap.onlineExists(channel)) continue;
      CRVROC onlineChannel  = crvChannelMap.online(channel);
      _ROCs.insert(onlineChannel.ROC());
    }
  }

  _totalEvents++;

  std::unique_ptr<artdaq::Fragments> fragments = std::make_unique<artdaq::Fragments>();
  auto crvDigiHandle = event.getHandle<mu2e::CrvDigiCollection>(_crvDigiTag);

  if(!crvDigiHandle.isValid())
  {
    if(_diagLevel > 0)
    {
      std::cout << "[CrvDigisToFragments::produce] Missing input CRV collection" << std::endl;
    }
    event.put(std::move(fragments));
    return;
  }

  auto const& crvDigis = *crvDigiHandle;
  _totalDigis += crvDigis.size();

  std::map<uint8_t, DTCLib::DTC_Event> dtcEvents;
  buildDtcEventsFromDigis(event, crvDigis, dtcEvents);

  for(auto& [dtcID, dtcEvent] : dtcEvents)
  {
    auto const eventBytes = dtcEvent.GetEventByteCount();

    if(_diagLevel > 1)
    {
      std::cout << "[CrvDigisToFragments::produce] DTC " << static_cast<int>(dtcID)
                << " eventBytes: " << eventBytes
                << ", subEventCount: " << dtcEvent.GetSubEventCount() << std::endl;
      for(size_t i = 0; i < dtcEvent.GetSubEvents().size(); ++i)
      {
        auto const& subEvt = dtcEvent.GetSubEvents()[i];
        std::cout << "  SubEvent " << i << " blockCount: " << subEvt.GetDataBlocks().size() << std::endl;
        for (size_t j = 0; j < subEvt.GetDataBlocks().size(); ++j)
        {
          auto const& block = subEvt.GetDataBlocks()[j];
          std::cout << "    Block " << j << " byteSize: " << block.byteSize << std::endl;
        }
      }
    }

    if(eventBytes == 0 || dtcEvent.GetSubEventCount() == 0) continue;

    std::vector<uint8_t> packed;
    packed.reserve(eventBytes);

    auto const* evHdr = dtcEvent.GetHeader();
    packed.insert(packed.end(), reinterpret_cast<uint8_t const*>(evHdr),
                  reinterpret_cast<uint8_t const*>(evHdr) + kDtcEventHeaderBytes);

    for(auto const& subEvt : dtcEvent.GetSubEvents())
    {
      auto const* subHdr = subEvt.GetHeader();
      packed.insert(packed.end(), reinterpret_cast<uint8_t const*>(subHdr),
                    reinterpret_cast<uint8_t const*>(subHdr) + kDtcSubEventHeaderBytes);

      for(auto const& block : subEvt.GetDataBlocks())
      {
        auto const* blk = reinterpret_cast<uint8_t const*>(block.GetRawBufferPointer());
        packed.insert(packed.end(), blk, blk + block.byteSize);
      }
    }

    if(_diagLevel > 0 && packed.size() != eventBytes)
    {
      std::cout << "[CrvDigisToFragments::produce] WARNING: DTC " << static_cast<int>(dtcID)
                << " packed size " << packed.size()
                << " differs from event header size " << eventBytes << std::endl;
    }

    auto fragPtr = artdaq::Fragment::FragmentBytes(packed.size());
    fragPtr->setUserType(mu2e::FragmentType::DTCEVT);
    fragPtr->setSequenceID(static_cast<artdaq::Fragment::sequence_id_t>(event.event()));
    fragPtr->setFragmentID(static_cast<artdaq::Fragment::fragment_id_t>(dtcID + _fragmentIdOffset));
    fragPtr->setTimestamp(static_cast<artdaq::Fragment::timestamp_t>(event.event()));
    if (!packed.empty()) std::memcpy(fragPtr->dataBeginBytes(), packed.data(), packed.size());

    fragments->emplace_back(std::move(*fragPtr));
  }

  if(_diagLevel > 0)
  {
    std::cout << "[CrvDigisToFragments::produce] Run " << event.run() << ", subrun "
              << event.subRun() << ", event " << event.event() << " has " << crvDigis.size()
              << " CrvDigis and created " << fragments->size() << " DTCEVT fragment(s)"
              << std::endl;
  }

  event.put(std::move(fragments));
}

void CrvDigisToFragments::endJob()
{
  if(_diagLevel > 0)
  {
    std::cout << "\n ----- [CrvDigisToFragments] Summary ----- " << std::endl;
    std::cout << "Total events: " << _totalEvents << std::endl;
    std::cout << "Total CrvDigis processed: " << _totalDigis << std::endl;
  }
}

} //namespace mu2e

DEFINE_ART_MODULE(mu2e::CrvDigisToFragments)
