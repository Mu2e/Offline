// ======================================================================
//
// Make CRVDigis from CRVDataDecoders
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/CRVConditions/inc/CRVOrdinal.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvDAQerror.hh"

#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include "artdaq-core-mu2e/Overlays/DTCEventFragment.hh"
#include "artdaq-core-mu2e/Overlays/Decoders/CRVDataDecoder.hh"
#include "artdaq-core/Data/Fragment.hh"

#include <iostream>
#include <string>
#include <memory>

namespace mu2e
{

class CrvDigisFromArtdaqFragmentsFEBII : public art::EDProducer
{
  public:
  struct Config
  {
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic Level")};
    //for currently wrongly encoded subevent headers
    fhicl::Atom<bool> useSubsystem0{fhicl::Name("useSubsystem0"), fhicl::Comment("consider subevents encoded with subsystem 0")};
  };

  explicit CrvDigisFromArtdaqFragmentsFEBII(const art::EDProducer::Table<Config>& config);
  ~CrvDigisFromArtdaqFragmentsFEBII() override {}
  void produce(art::Event&) override;

  private:
  int                                      _diagLevel;
  bool                                     _useSubsystem0;
  mu2e::ProditionsHandle<mu2e::CRVOrdinal> _channelMap_h;
};

CrvDigisFromArtdaqFragmentsFEBII::CrvDigisFromArtdaqFragmentsFEBII(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}, _diagLevel(config().diagLevel()), _useSubsystem0(config().useSubsystem0())
{
  produces<mu2e::CrvDigiCollection>();
  produces<mu2e::CrvDigiCollection>("NZS");
  produces<mu2e::CrvDAQerrorCollection>();
}

void CrvDigisFromArtdaqFragmentsFEBII::produce(art::Event& event)
{
  // Collection of CrvDigis for the event
  std::unique_ptr<mu2e::CrvDigiCollection> crvDigis(new mu2e::CrvDigiCollection);
  std::unique_ptr<mu2e::CrvDigiCollection> crvDigisNZS(new mu2e::CrvDigiCollection);
  std::unique_ptr<mu2e::CrvDAQerrorCollection> crvDaqErrors(new mu2e::CrvDAQerrorCollection);

  auto const& channelMap = _channelMap_h.get(event.id());

  art::EventNumber_t eventNumber = event.event();
  if(_diagLevel>1)
  {
    std::cout << std::dec << "Run/Subrun/Event: " << event.run() << "/" << event.subRun() << "/" << eventNumber << std::endl;
  }

  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;
  fragmentHandles = event.getMany<std::vector<artdaq::Fragment>>();

  size_t totalSize = 0;
  size_t iFragment = 0;
  for(const auto& fragmentHandle : fragmentHandles)
  {
    if(!fragmentHandle.isValid()) continue;

    for(auto artFragment : *fragmentHandle)
    {
      if(artFragment.type()!=mu2e::FragmentType::DTCEVT) continue;

      mu2e::DTCEventFragment dtcEventFragment(artFragment);
      ++iFragment;

      //collect sub events
      auto dtcSubEvents = dtcEventFragment.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_CRV);
      if(_useSubsystem0)
      {
        auto dtcSubEventsTmp = dtcEventFragment.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker);  //currently wrongly encoded in the DTC Subevent header
        dtcSubEvents.insert(dtcSubEvents.end(),dtcSubEventsTmp.begin(),dtcSubEventsTmp.end()); //temporarily add fragments with wrongly encoded headers
      }

      //check for errors
      const DTCLib::DTC_Event &dtcEvent =  dtcEventFragment.getData();
      size_t expectedSize = dtcEvent.GetEventByteCount();
      size_t actualSize = sizeof(DTCLib::DTC_EventHeader);
      for(size_t iSubEvent=0; iSubEvent<dtcSubEvents.size(); ++iSubEvent) actualSize+=dtcSubEvents.at(iSubEvent).GetSubEventByteCount();
      if(_diagLevel>1)
      {
        std::cout << "Fragment index: " << iFragment << "      expected event size: " << expectedSize << ", actual event size: " << actualSize << std::endl;
      }
      if(expectedSize!=actualSize)
      {
        std::cerr << "mismatch between expected event size and actual event size!" << std::endl;
        crvDaqErrors->emplace_back(mu2e::CrvDAQerrorCode::byteCountMismatch,iFragment,0,0,0);
      }

      //loop through subevents
      for(size_t iSubEvent=0; auto& dtcSubEvent : dtcSubEvents)
      {
        mu2e::CRVDataDecoder decoder(dtcSubEvent);
        for(size_t iDataBlock = 0; iDataBlock < decoder.block_count(); ++iDataBlock)
        {
          totalSize+=decoder.blockSizeBytes(iDataBlock);

          auto block = decoder.dataAtBlockIndex(iDataBlock);
          if(block == nullptr)
          {
            std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
            std::cerr << "Unable to retrieve data block." << std::endl;
            crvDaqErrors->emplace_back(mu2e::CrvDAQerrorCode::unableToGetDataBlock,iFragment,iSubEvent,iDataBlock,0);
            continue;
          }

          auto header = block->GetHeader();
          if(!header->isValid())
          {
            std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
            std::cerr << "CRV packet is not valid." << std::endl;
            std::cerr << "sub system ID: "<<(uint16_t)header->GetSubsystemID()<<" packet count: "<<header->GetPacketCount() << std::endl;
            crvDaqErrors->emplace_back(mu2e::CrvDAQerrorCode::invalidPacket,iFragment,iSubEvent,iDataBlock,header->GetPacketCount());
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
            crvDaqErrors->emplace_back(mu2e::CrvDAQerrorCode::wrongSubsystemID,iFragment,iSubEvent,iDataBlock,header->GetPacketCount());
            continue;
          }

          if(_diagLevel>1)
          {
            std::cout << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
            std::cout << "packet count: " << header->GetPacketCount() << std::endl;
          }
          if(header->GetPacketCount()>0)
          {
            auto crvRocHeader = decoder.GetCRVROCStatusPacket(iDataBlock);
            if(crvRocHeader == nullptr)
            {
              std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
              std::cerr << "Error retrieving CRV ROC Status Packet" << std::endl;
              crvDaqErrors->emplace_back(mu2e::CrvDAQerrorCode::errorUnpackingStatusPacket,iFragment,iSubEvent,iDataBlock,header->GetPacketCount());
              continue;
            }

            std::vector<mu2e::CRVDataDecoder::CRVHitFEBII> crvHits;
            if(!decoder.GetCRVHitsFEBII(iDataBlock, crvHits))
            {
              std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
              std::cerr << "Error unpacking of CRV Hits" << std::endl;
              crvDaqErrors->emplace_back(mu2e::CrvDAQerrorCode::errorUnpackingCrvHits,iFragment,iSubEvent,iDataBlock,header->GetPacketCount());
              break;
            }
            for(auto const& crvHit : crvHits)
            {
              const auto& crvHitInfo = crvHit.first;
              const auto& waveform = crvHit.second;

              uint16_t rocID = crvRocHeader->ControllerID+1; // FIXME ROC IDs between 1 and 17  //also header->GetLinkID()
              uint16_t rocPort = crvHitInfo.portNumber+1; //FIXME Port numbers beween 1 and 24
              uint16_t febChannel = crvHitInfo.febChannel;
              mu2e::CRVROC onlineChannel(rocID, rocPort, febChannel);

              uint16_t offlineChannel = channelMap.offline(onlineChannel);
              int crvBarIndex = offlineChannel / 4;
              int SiPMNumber = offlineChannel % 4;

              std::vector<int16_t> adc;
              adc.insert(adc.end(),waveform.begin(),waveform.end());  //convert from uint16_t to int16_t
              for(size_t i=0; i<waveform.size(); ++i) {if((adc[i] & 0x800) == 0x800) adc[i]=(int16_t)(adc[i] | 0xF000);}  //to handle negative numbers stored in 12bit ADC samples
              crvDigis->emplace_back(adc, crvHitInfo.hitTime, false, mu2e::CRSScintillatorBarIndex(crvBarIndex), SiPMNumber);
              crvDigisNZS->emplace_back(adc, crvHitInfo.hitTime, true, mu2e::CRSScintillatorBarIndex(crvBarIndex), SiPMNumber);  //temporary solution until we get the FEB-II
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

                uint16_t rocID = crvRocHeader->ControllerID + 1; // FIXME  //ROC IDs are between 1 and 17
                uint16_t rocPort = crvHitInfo.portNumber+1; //FIXME Port numbers beween 1 and 24
                uint16_t febChannel = crvHitInfo.febChannel;
                mu2e::CRVROC onlineChannel(rocID, rocPort, febChannel);

                uint16_t offlineChannel = channelMap.offline(onlineChannel);
                int crvBarIndex = offlineChannel / 4;
                int SiPMNumber = offlineChannel % 4;

                std::cout << "ROCID (increased by 1 to match the Online/Offline-Channel Map) " << rocID
                          << "   rocPort " << rocPort << "   febChannel " << febChannel
                          << "   crvBarIndex " << crvBarIndex << "   SiPMNumber " << SiPMNumber
                          << std::endl;
                std::cout << "TDC: " << crvHitInfo.hitTime << std::endl;
                std::cout << "nSamples " << waveform.size() << "  ";
                std::cout << "Waveform: {";
                for(size_t i = 0; i < waveform.size(); ++i) std::cout << "  " << waveform.at(i);
                std::cout << "}" << std::endl;
                std::cout << std::endl;
              } // loop over hits
            }   // debug output
          }     // packet count > 0
        }       // loop over data blocks
      }         // loop over subEvents(=decoders)
    }           // loop over fragments
  }             // loop over fragment handles

  if(_diagLevel>1) std::cout << "Total Size: " << totalSize << std::endl;

  // Store the crv digis in the event
  event.put(std::move(crvDigis));
  event.put(std::move(crvDigisNZS),"NZS");
  event.put(std::move(crvDaqErrors));
}

} //namespace mu2e

using mu2e::CrvDigisFromArtdaqFragmentsFEBII;
DEFINE_ART_MODULE(CrvDigisFromArtdaqFragmentsFEBII)
