// ======================================================================
//
// Make CRVDigis from CRVDataDecoders
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/CRVConditions/inc/CRVOrdinal.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvDAQerror.hh"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Data/CRVDataDecoder.hh"
#include <artdaq-core/Data/Fragment.hh>

#include <iostream>
#include <string>
#include <memory>

namespace art
{
  class CRVDigisFromDataDecoders;
}

using art::CRVDigisFromDataDecoders;

class art::CRVDigisFromDataDecoders : public EDProducer
{
  public:
  struct Config
  {
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic Level")};
    fhicl::Atom<art::InputTag> CRVDataDecodersTag{fhicl::Name("crvTag"), fhicl::Comment("crv Fragments Tag")};
  };

  explicit CRVDigisFromDataDecoders(const art::EDProducer::Table<Config>& config);
  ~CRVDigisFromDataDecoders() override {}

  void produce(art::Event&) override;

  private:
  int                                      _diagLevel;
  art::InputTag                            _CRVDataDecodersTag;
  mu2e::ProditionsHandle<mu2e::CRVOrdinal> _channelMap_h;

};

CRVDigisFromDataDecoders::CRVDigisFromDataDecoders(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}, _diagLevel(config().diagLevel()), _CRVDataDecodersTag(config().CRVDataDecodersTag())
{
  produces<mu2e::CrvDigiCollection>();
  produces<mu2e::CrvDigiCollection>("NZS");
  produces<mu2e::CrvDAQerrorCollection>();
}

void CRVDigisFromDataDecoders::produce(Event& event)
{
  art::EventNumber_t eventNumber = event.event();

  auto CRVDataDecoders = event.getValidHandle<std::vector<mu2e::CRVDataDecoder> >(_CRVDataDecodersTag);
  size_t nSubEvents = CRVDataDecoders->size();

  if(_diagLevel>1)
  {
    std::cout << std::dec << "Run/Subrun/Event: " << event.run() << "/" << event.subRun() << "/" << eventNumber << std::endl;

    size_t totalSize = 0;
    for(const auto &frag : *CRVDataDecoders)
    {
      for(size_t i = 0; i < frag.block_count(); ++i)
      {
        auto size = frag.blockSizeBytes(i);
        totalSize += size;
      }
    }

    std::cout << "#SubEvents: " << nSubEvents << std::endl;
    std::cout << "Total Size: " << totalSize << " bytes." << std::endl;
  }

  // Collection of CrvDigis for the event
  std::unique_ptr<mu2e::CrvDigiCollection> crv_digis(new mu2e::CrvDigiCollection);
  std::unique_ptr<mu2e::CrvDigiCollection> crv_digis_NZS(new mu2e::CrvDigiCollection);
  std::unique_ptr<mu2e::CrvDAQerrorCollection> crv_daq_errors(new mu2e::CrvDAQerrorCollection);
  auto const& channelMap = _channelMap_h.get(event.id());

  // Loop over the CRV fragments
  for(size_t iSubEvent = 0; iSubEvent < nSubEvents; ++iSubEvent)
  {
    const mu2e::CRVDataDecoder& CRVDataDecoder((*CRVDataDecoders)[iSubEvent]);
    CRVDataDecoder.setup_event();

    for(size_t iDataBlock = 0; iDataBlock < CRVDataDecoder.block_count(); ++iDataBlock)
    {
      auto block = CRVDataDecoder.dataAtBlockIndex(iDataBlock);
      if(block == nullptr)
      {
        std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
        std::cerr << "Unable to retrieve data block." << std::endl;
        crv_daq_errors->emplace_back(mu2e::CrvDAQerrorCode::unableToGetDataBlock,iSubEvent,iDataBlock,0);
        continue;
      }
      auto header = block->GetHeader();

      if(!header->isValid())
      {
        std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
        std::cerr << "CRV packet is not valid." << std::endl;
        std::cerr << "sub system ID: "<<(uint16_t)header->GetSubsystemID()<<" packet count: "<<header->GetPacketCount() << std::endl;
        crv_daq_errors->emplace_back(mu2e::CrvDAQerrorCode::invalidPacket,iSubEvent,iDataBlock,header->GetPacketCount());
        continue;
      }

      if(header->GetSubsystemID() != DTCLib::DTC_Subsystem::DTC_Subsystem_CRV)
      {
        if(_diagLevel>0)
        {
          std::cout << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
          std::cout << "CRV packet does not have subsystem ID 2." << std::endl;
          std::cout << "sub system ID: "<<(uint16_t)header->GetSubsystemID()<<" packet count: "<<header->GetPacketCount() << std::endl;
        }
        crv_daq_errors->emplace_back(mu2e::CrvDAQerrorCode::wrongSubsystemID,iSubEvent,iDataBlock,header->GetPacketCount());
        continue;
      }

      if(_diagLevel>1)
      {
        std::cout << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
        std::cout << "packet count: " << header->GetPacketCount() << std::endl;
      }
      if(header->GetPacketCount() > 0)
      {
        auto crvRocHeader = CRVDataDecoder.GetCRVROCStatusPacket(iDataBlock);
        if(crvRocHeader == nullptr)
        {
          std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
          std::cerr << "Error retrieving CRV ROC Status Packet" << std::endl;
          crv_daq_errors->emplace_back(mu2e::CrvDAQerrorCode::errorUnpackingStatusPacket,iSubEvent,iDataBlock,header->GetPacketCount());
          continue;
        }

        std::vector<mu2e::CRVDataDecoder::CRVHit> crvHits;
        if(!CRVDataDecoder.GetCRVHits(iDataBlock, crvHits))
        {
          std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
          std::cerr << "Error unpacking of CRV Hits" << std::endl;
          crv_daq_errors->emplace_back(mu2e::CrvDAQerrorCode::errorUnpackingCrvHits,iSubEvent,iDataBlock,header->GetPacketCount());
          break;
        }
        for(auto const& crvHit : crvHits)
        {
          const auto& crvHitInfo = crvHit.first;
          const auto& waveform = crvHit.second;

          uint16_t rocID = crvHitInfo.controllerNumber + 1; // FIXME ROC IDs between 1 and 17   //also: header->GetLinkID()+1;
          uint16_t rocPort = crvHitInfo.portNumber;
          uint16_t febChannel = crvHitInfo.febChannel;
          mu2e::CRVROC onlineChannel(rocID, rocPort, febChannel);

          uint16_t offlineChannel = channelMap.offline(onlineChannel);
          int crvBarIndex = offlineChannel / 4;
          int SiPMNumber = offlineChannel % 4;

          std::vector<int16_t> adc;
          adc.resize(waveform.size());
          for(size_t i=0; i<waveform.size(); ++i) adc[i]=waveform.at(i).ADC;
          for(size_t i=0; i<waveform.size(); ++i) {if((adc[i] & 0x800) == 0x800) adc[i]=(int16_t)(adc[i] | 0xF000);}  //to handle negative numbers stored in 12bit ADC samples
          crv_digis->emplace_back(adc, crvHitInfo.HitTime, false, mu2e::CRSScintillatorBarIndex(crvBarIndex), SiPMNumber);
          crv_digis_NZS->emplace_back(adc, crvHitInfo.HitTime, true, mu2e::CRSScintillatorBarIndex(crvBarIndex), SiPMNumber);  //temporary solution until we get the FEB-II
        } // loop over all crvHits

        if(_diagLevel>2)
        {
          std::cout << "EventWindowTag (TDC header): " << header->GetEventWindowTag().GetEventWindowTag(true) << std::endl;
          std::cout << "SubsystemID: " << (uint16_t)header->GetSubsystemID() << std::endl;
          std::cout << "DTCID: " << (uint16_t)header->GetID() << std::endl;
          std::cout << "ROCID (TDC header): " << (uint16_t)header->GetLinkID() << std::endl;
          std::cout << "packetCount: " << header->GetPacketCount() << std::endl;
          std::cout << "EVB mode: " << (uint16_t)header->GetEVBMode() << std::endl;
          std::cout << "TriggerCount: " << crvRocHeader->TriggerCount << std::endl;
          std::cout << "ActiveFEBFlags: " << crvRocHeader->GetActiveFEBFlags() << std::endl;
          std::cout << "MicroBunchStatus: " << std::hex << (uint16_t)crvRocHeader->MicroBunchStatus << std::dec << std::endl;
          std::cout << "ROCID (ROC header): " << (uint16_t)crvRocHeader->ControllerID  << std::endl;
          std::cout << "EventWindowTag (ROC header): " << crvRocHeader->GetEventWindowTag() << std::endl;
        }

        if(_diagLevel>3)
        {
          for(auto const& crvHit : crvHits)
          {
            const auto& crvHitInfo = crvHit.first;
            const auto& waveform = crvHit.second;

            uint16_t rocID = crvHitInfo.controllerNumber + 1; // FIXME  //ROC IDs are between 1 and 17
            uint16_t rocPort = crvHitInfo.portNumber;
            uint16_t febChannel = crvHitInfo.febChannel;
            mu2e::CRVROC onlineChannel(rocID, rocPort, febChannel);

            uint16_t offlineChannel = channelMap.offline(onlineChannel);
            int crvBarIndex = offlineChannel / 4;
            int SiPMNumber = offlineChannel % 4;

            std::cout << "ROCID (increased by 1 to match the Online/Offline-Channel Map) " << rocID
                      << "   rocPort " << rocPort << "   febChannel " << febChannel
                      << "   crvBarIndex " << crvBarIndex << "   SiPMNumber " << SiPMNumber
                      << std::endl;
            std::cout << "TDC: " << crvHitInfo.HitTime << std::endl;
            std::cout << "nSamples " << crvHitInfo.NumSamples << "  ";
            std::cout << "Waveform: {";
            for(size_t i = 0; i < crvHitInfo.NumSamples; i++)
              std::cout << "  " << waveform.at(i).ADC;
            std::cout << "}" << std::endl;
            std::cout << std::endl;
          } // loop over hits
        }   // debug output
      }     // end parsing CRV DataBlocks
    }       // loop over DataBlocks within CRVDataDecoders
  }         // Close loop over fragments

  // Store the crv digis in the event
  event.put(std::move(crv_digis));
  event.put(std::move(crv_digis_NZS),"NZS");
  event.put(std::move(crv_daq_errors));
}

DEFINE_ART_MODULE(CRVDigisFromDataDecoders)
